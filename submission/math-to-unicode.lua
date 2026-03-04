-- =============================================================================
-- FILTER: math-to-unicode.lua
-- PURPOSE: Pandoc Lua filter to convert LaTeX math to Unicode text for .docx
--          output. Maps Greek letters, math operators, superscripts, subscripts,
--          and common symbols to their Unicode equivalents.
-- USAGE:   pandoc --lua-filter=math-to-unicode.lua input.tex -o output.docx
-- PROJECT: Core-PAM (Memorial CorePAM)
-- =============================================================================

-- Greek letter mapping (lowercase + uppercase)
local greek = {
  ["\\alpha"]   = "\u{03B1}",
  ["\\beta"]    = "\u{03B2}",
  ["\\gamma"]   = "\u{03B3}",
  ["\\delta"]   = "\u{03B4}",
  ["\\epsilon"] = "\u{03B5}",
  ["\\varepsilon"] = "\u{03B5}",
  ["\\zeta"]    = "\u{03B6}",
  ["\\eta"]     = "\u{03B7}",
  ["\\theta"]   = "\u{03B8}",
  ["\\iota"]    = "\u{03B9}",
  ["\\kappa"]   = "\u{03BA}",
  ["\\lambda"]  = "\u{03BB}",
  ["\\mu"]      = "\u{03BC}",
  ["\\nu"]      = "\u{03BD}",
  ["\\xi"]      = "\u{03BE}",
  ["\\pi"]      = "\u{03C0}",
  ["\\rho"]     = "\u{03C1}",
  ["\\sigma"]   = "\u{03C3}",
  ["\\tau"]     = "\u{03C4}",
  ["\\upsilon"] = "\u{03C5}",
  ["\\phi"]     = "\u{03C6}",
  ["\\varphi"]  = "\u{03C6}",
  ["\\chi"]     = "\u{03C7}",
  ["\\psi"]     = "\u{03C8}",
  ["\\omega"]   = "\u{03C9}",
  ["\\Alpha"]   = "\u{0391}",
  ["\\Beta"]    = "\u{0392}",
  ["\\Gamma"]   = "\u{0393}",
  ["\\Delta"]   = "\u{0394}",
  ["\\Theta"]   = "\u{0398}",
  ["\\Lambda"]  = "\u{039B}",
  ["\\Pi"]      = "\u{03A0}",
  ["\\Sigma"]   = "\u{03A3}",
  ["\\Phi"]     = "\u{03A6}",
  ["\\Psi"]     = "\u{03A8}",
  ["\\Omega"]   = "\u{03A9}",
}

-- Superscript digits and common chars
local superscripts = {
  ["0"] = "\u{2070}", ["1"] = "\u{00B9}", ["2"] = "\u{00B2}",
  ["3"] = "\u{00B3}", ["4"] = "\u{2074}", ["5"] = "\u{2075}",
  ["6"] = "\u{2076}", ["7"] = "\u{2077}", ["8"] = "\u{2078}",
  ["9"] = "\u{2079}", ["+"] = "\u{207A}", ["-"] = "\u{207B}",
  ["="] = "\u{207C}", ["("] = "\u{207D}", [")"] = "\u{207E}",
  ["n"] = "\u{207F}", ["i"] = "\u{2071}",
}

-- Subscript digits
local subscripts = {
  ["0"] = "\u{2080}", ["1"] = "\u{2081}", ["2"] = "\u{2082}",
  ["3"] = "\u{2083}", ["4"] = "\u{2084}", ["5"] = "\u{2085}",
  ["6"] = "\u{2086}", ["7"] = "\u{2087}", ["8"] = "\u{2088}",
  ["9"] = "\u{2089}", ["+"] = "\u{208A}", ["-"] = "\u{208B}",
  ["="] = "\u{208C}", ["("] = "\u{208D}", [")"] = "\u{208E}",
  ["i"] = "\u{1D62}", ["j"] = "\u{2C7C}",
}

-- Math operators and symbols
local symbols = {
  ["\\times"]    = "\u{00D7}",
  ["\\cdot"]     = "\u{00B7}",
  ["\\pm"]       = "\u{00B1}",
  ["\\mp"]       = "\u{2213}",
  ["\\leq"]      = "\u{2264}",
  ["\\le"]       = "\u{2264}",
  ["\\geq"]      = "\u{2265}",
  ["\\ge"]       = "\u{2265}",
  ["\\neq"]      = "\u{2260}",
  ["\\ne"]       = "\u{2260}",
  ["\\approx"]   = "\u{2248}",
  ["\\sim"]      = "\u{223C}",
  ["\\infty"]    = "\u{221E}",
  ["\\sum"]      = "\u{2211}",
  ["\\prod"]     = "\u{220F}",
  ["\\int"]      = "\u{222B}",
  ["\\partial"]  = "\u{2202}",
  ["\\nabla"]    = "\u{2207}",
  ["\\in"]       = "\u{2208}",
  ["\\notin"]    = "\u{2209}",
  ["\\subset"]   = "\u{2282}",
  ["\\supset"]   = "\u{2283}",
  ["\\cup"]      = "\u{222A}",
  ["\\cap"]      = "\u{2229}",
  ["\\forall"]   = "\u{2200}",
  ["\\exists"]   = "\u{2203}",
  ["\\rightarrow"] = "\u{2192}",
  ["\\leftarrow"]  = "\u{2190}",
  ["\\Rightarrow"] = "\u{21D2}",
  ["\\Leftarrow"]  = "\u{21D0}",
  ["\\ldots"]    = "\u{2026}",
  ["\\cdots"]    = "\u{22EF}",
  ["\\quad"]     = " ",
  ["\\qquad"]    = "  ",
  ["\\,"]        = "\u{2009}",  -- thin space
  ["\\;"]        = "\u{2005}",  -- medium space
  ["\\!"]        = "",           -- negative thin space (remove)
  ["\\text"]     = "",           -- will be handled specially
}

-- Convert a string of characters to superscript Unicode
local function to_superscript(s)
  local result = ""
  for c in s:gmatch(".") do
    result = result .. (superscripts[c] or c)
  end
  return result
end

-- Convert a string of characters to subscript Unicode
local function to_subscript(s)
  local result = ""
  for c in s:gmatch(".") do
    result = result .. (subscripts[c] or c)
  end
  return result
end

-- Main math-to-unicode converter
local function math_to_unicode(s)
  -- Remove \mathrm{}, \text{}, \textit{} wrappers (keep content)
  s = s:gsub("\\mathrm%s*{([^}]*)}", "%1")
  s = s:gsub("\\text%s*{([^}]*)}", "%1")
  s = s:gsub("\\textit%s*{([^}]*)}", "%1")
  s = s:gsub("\\mathbf%s*{([^}]*)}", "%1")
  s = s:gsub("\\boldsymbol%s*{([^}]*)}", "%1")

  -- Handle \frac{a}{b} -> a/b
  s = s:gsub("\\frac%s*{([^}]*)}{([^}]*)}", "%1/%2")

  -- Handle superscripts: ^{...} and ^x
  s = s:gsub("%^{([^}]*)}", function(content)
    return to_superscript(content)
  end)
  s = s:gsub("%^(%w)", function(c)
    return to_superscript(c)
  end)

  -- Handle subscripts: _{...} and _x
  s = s:gsub("_{([^}]*)}", function(content)
    return to_subscript(content)
  end)
  s = s:gsub("_(%w)", function(c)
    return to_subscript(c)
  end)

  -- Replace symbols (longer patterns first to avoid partial matches)
  local sorted_keys = {}
  for k, _ in pairs(symbols) do sorted_keys[#sorted_keys + 1] = k end
  table.sort(sorted_keys, function(a, b) return #a > #b end)
  for _, k in ipairs(sorted_keys) do
    s = s:gsub(k:gsub("\\", "\\\\"), symbols[k])
  end

  -- Replace Greek letters (longer patterns first)
  local greek_keys = {}
  for k, _ in pairs(greek) do greek_keys[#greek_keys + 1] = k end
  table.sort(greek_keys, function(a, b) return #a > #b end)
  for _, k in ipairs(greek_keys) do
    s = s:gsub(k:gsub("\\", "\\\\"), greek[k])
  end

  -- Remove remaining braces
  s = s:gsub("[{}]", "")

  -- Clean up multiple spaces
  s = s:gsub("  +", " ")

  return s
end

-- Pandoc filter: convert Math elements to Str
function Math(el)
  local converted = math_to_unicode(el.text)
  return pandoc.Str(converted)
end
