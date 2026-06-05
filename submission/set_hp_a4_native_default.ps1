$ErrorActionPreference = 'Stop'

$printerName = 'HP Color LaserJet Pro MFP 4303'

Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;

public static class NativePrinterDefaults
{
    [StructLayout(LayoutKind.Sequential, CharSet = CharSet.Unicode)]
    public struct DEVMODE
    {
        [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 32)]
        public string dmDeviceName;
        public short dmSpecVersion;
        public short dmDriverVersion;
        public short dmSize;
        public short dmDriverExtra;
        public int dmFields;
        public short dmOrientation;
        public short dmPaperSize;
        public short dmPaperLength;
        public short dmPaperWidth;
        public short dmScale;
        public short dmCopies;
        public short dmDefaultSource;
        public short dmPrintQuality;
        public short dmColor;
        public short dmDuplex;
        public short dmYResolution;
        public short dmTTOption;
        public short dmCollate;
        [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 32)]
        public string dmFormName;
        public short dmLogPixels;
        public int dmBitsPerPel;
        public int dmPelsWidth;
        public int dmPelsHeight;
        public int dmDisplayFlags;
        public int dmDisplayFrequency;
        public int dmICMMethod;
        public int dmICMIntent;
        public int dmMediaType;
        public int dmDitherType;
        public int dmReserved1;
        public int dmReserved2;
        public int dmPanningWidth;
        public int dmPanningHeight;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct PRINTER_INFO_9
    {
        public IntPtr pDevMode;
    }

    [DllImport("winspool.drv", SetLastError = true, CharSet = CharSet.Unicode)]
    public static extern bool OpenPrinter(string pPrinterName, out IntPtr phPrinter, IntPtr pDefault);

    [DllImport("winspool.drv", SetLastError = true)]
    public static extern bool ClosePrinter(IntPtr hPrinter);

    [DllImport("winspool.drv", SetLastError = true, CharSet = CharSet.Unicode)]
    public static extern int DocumentProperties(
        IntPtr hwnd,
        IntPtr hPrinter,
        string pDeviceName,
        IntPtr pDevModeOutput,
        IntPtr pDevModeInput,
        int fMode);

    [DllImport("winspool.drv", SetLastError = true)]
    public static extern bool SetPrinter(IntPtr hPrinter, int Level, IntPtr pPrinter, int Command);
}
'@

$DM_OUT_BUFFER = 0x00000002
$DM_IN_BUFFER = 0x00000008

$DM_ORIENTATION = 0x00000001
$DM_PAPERSIZE = 0x00000002
$DM_COLOR = 0x00000800
$DM_DUPLEX = 0x00001000
$DM_COLLATE = 0x00008000

$DMORIENT_PORTRAIT = 1
$DMPAPER_A4 = 9
$DMCOLOR_COLOR = 2
$DMDUP_VERTICAL_LONG_EDGE = 2
$DMCOLLATE_TRUE = 1

$hPrinter = [IntPtr]::Zero
if (-not [NativePrinterDefaults]::OpenPrinter($printerName, [ref]$hPrinter, [IntPtr]::Zero)) {
    throw "OpenPrinter failed: $([Runtime.InteropServices.Marshal]::GetLastWin32Error())"
}

try {
    $needed = [NativePrinterDefaults]::DocumentProperties([IntPtr]::Zero, $hPrinter, $printerName, [IntPtr]::Zero, [IntPtr]::Zero, 0)
    if ($needed -le 0) {
        throw "DocumentProperties(size) failed: $([Runtime.InteropServices.Marshal]::GetLastWin32Error())"
    }

    $devModePtr = [Runtime.InteropServices.Marshal]::AllocHGlobal($needed)
    try {
        $result = [NativePrinterDefaults]::DocumentProperties([IntPtr]::Zero, $hPrinter, $printerName, $devModePtr, [IntPtr]::Zero, $DM_OUT_BUFFER)
        if ($result -lt 0) {
            throw "DocumentProperties(read) failed: $([Runtime.InteropServices.Marshal]::GetLastWin32Error())"
        }

        $devMode = [Runtime.InteropServices.Marshal]::PtrToStructure($devModePtr, [type][NativePrinterDefaults+DEVMODE])
        $before = [pscustomobject]@{
            PaperSize = $devMode.dmPaperSize
            Duplex = $devMode.dmDuplex
            Color = $devMode.dmColor
            Collate = $devMode.dmCollate
            Orientation = $devMode.dmOrientation
            Fields = ('0x{0:X}' -f $devMode.dmFields)
        }

        $devMode.dmFields = $devMode.dmFields -bor $DM_ORIENTATION -bor $DM_PAPERSIZE -bor $DM_COLOR -bor $DM_DUPLEX -bor $DM_COLLATE
        $devMode.dmOrientation = $DMORIENT_PORTRAIT
        $devMode.dmPaperSize = $DMPAPER_A4
        $devMode.dmColor = $DMCOLOR_COLOR
        $devMode.dmDuplex = $DMDUP_VERTICAL_LONG_EDGE
        $devMode.dmCollate = $DMCOLLATE_TRUE
        [Runtime.InteropServices.Marshal]::StructureToPtr($devMode, $devModePtr, $false)

        $result = [NativePrinterDefaults]::DocumentProperties([IntPtr]::Zero, $hPrinter, $printerName, $devModePtr, $devModePtr, ($DM_IN_BUFFER -bor $DM_OUT_BUFFER))
        if ($result -lt 0) {
            throw "DocumentProperties(validate) failed: $([Runtime.InteropServices.Marshal]::GetLastWin32Error())"
        }

        $info = New-Object NativePrinterDefaults+PRINTER_INFO_9
        $info.pDevMode = $devModePtr
        $infoPtr = [Runtime.InteropServices.Marshal]::AllocHGlobal([Runtime.InteropServices.Marshal]::SizeOf($info))
        try {
            [Runtime.InteropServices.Marshal]::StructureToPtr($info, $infoPtr, $false)
            if (-not [NativePrinterDefaults]::SetPrinter($hPrinter, 9, $infoPtr, 0)) {
                throw "SetPrinter(level 9) failed: $([Runtime.InteropServices.Marshal]::GetLastWin32Error())"
            }
        }
        finally {
            [Runtime.InteropServices.Marshal]::FreeHGlobal($infoPtr)
        }

        $result = [NativePrinterDefaults]::DocumentProperties([IntPtr]::Zero, $hPrinter, $printerName, $devModePtr, [IntPtr]::Zero, $DM_OUT_BUFFER)
        if ($result -lt 0) {
            throw "DocumentProperties(re-read) failed: $([Runtime.InteropServices.Marshal]::GetLastWin32Error())"
        }
        $afterDevMode = [Runtime.InteropServices.Marshal]::PtrToStructure($devModePtr, [type][NativePrinterDefaults+DEVMODE])
        $after = [pscustomobject]@{
            PaperSize = $afterDevMode.dmPaperSize
            Duplex = $afterDevMode.dmDuplex
            Color = $afterDevMode.dmColor
            Collate = $afterDevMode.dmCollate
            Orientation = $afterDevMode.dmOrientation
            Fields = ('0x{0:X}' -f $afterDevMode.dmFields)
        }

        [pscustomobject]@{
            Printer = $printerName
            Before = $before
            After = $after
            ExpectedPaperSize = 'A4 = DMPAPER_A4 = 9'
            ExpectedDuplex = 'Long edge = DMDUP_VERTICAL = 2'
        }
    }
    finally {
        [Runtime.InteropServices.Marshal]::FreeHGlobal($devModePtr)
    }
}
finally {
    [void][NativePrinterDefaults]::ClosePrinter($hPrinter)
}

