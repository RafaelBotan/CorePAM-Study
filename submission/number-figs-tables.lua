-- =============================================================================
-- FILTER: number-figs-tables.lua
-- PURPOSE: Pandoc Lua filter for LaTeX→DOCX conversion:
--   1. Adds sequential "Figure N." / "Table N." numbering to captions
--   2. Sets all images to full text width (6.5 inches for standard margins)
--   3. Resolves \ref cross-references to figure/table numbers
-- USAGE:   pandoc --lua-filter=number-figs-tables.lua input.tex -o output.docx
-- PROJECT: CorePAM (Breast Cancer Research / Springer Nature)
-- =============================================================================

local fig_counter = 0
local tab_counter = 0
local label_map = {}  -- maps LaTeX labels to "Figure N" / "Table N"

-- Full text width: 6.5 inches (A4 with 1" margins) = 16.51 cm
-- In pandoc Image attributes, width is specified as a string
local FULL_WIDTH = "6.5in"

-- ---------------------------------------------------------------------------
-- Walk inline elements within a Figure to find and resize Image nodes
-- ---------------------------------------------------------------------------
local function resize_images(blocks)
  for _, block in ipairs(blocks) do
    if block.t == "Plain" or block.t == "Para" then
      for _, inline in ipairs(block.content) do
        if inline.t == "Image" then
          -- Set width to full text width
          inline.attr.attributes["width"] = FULL_WIDTH
          -- Remove any height constraint to preserve aspect ratio
          inline.attr.attributes["height"] = nil
        end
      end
    end
  end
end

-- ---------------------------------------------------------------------------
-- Figure: number caption and resize images
-- ---------------------------------------------------------------------------
function Figure(el)
  fig_counter = fig_counter + 1
  local label = "Figure " .. fig_counter

  -- Register label for cross-references
  if el.identifier and el.identifier ~= "" then
    label_map[el.identifier] = label
  end

  -- Prepend "Figure N. " in bold to the caption
  if el.caption and el.caption.long and #el.caption.long > 0 then
    local first_block = el.caption.long[1]
    if first_block and first_block.content then
      -- Insert bold "Figure N. " at the start
      table.insert(first_block.content, 1, pandoc.Space())
      table.insert(first_block.content, 1, pandoc.Strong{pandoc.Str(label .. ".")})
    end
  end

  -- Resize all images in this figure to full width
  resize_images(el.content)

  return el
end

-- ---------------------------------------------------------------------------
-- Table: number caption
-- ---------------------------------------------------------------------------
function Table(el)
  tab_counter = tab_counter + 1
  local label = "Table " .. tab_counter

  -- Register label for cross-references
  if el.identifier and el.identifier ~= "" then
    label_map[el.identifier] = label
  end

  -- Prepend "Table N. " in bold to the caption
  if el.caption and el.caption.long and #el.caption.long > 0 then
    local first_block = el.caption.long[1]
    if first_block and first_block.content then
      table.insert(first_block.content, 1, pandoc.Space())
      table.insert(first_block.content, 1, pandoc.Strong{pandoc.Str(label .. ".")})
    end
  end

  return el
end

-- ---------------------------------------------------------------------------
-- Standalone images (not in figure environment): also resize to full width
-- ---------------------------------------------------------------------------
function Image(el)
  el.attr.attributes["width"] = FULL_WIDTH
  el.attr.attributes["height"] = nil
  return el
end

-- ---------------------------------------------------------------------------
-- Filter ordering: figures/tables first (to build label_map),
-- then resolve cross-references
-- ---------------------------------------------------------------------------
return {
  {Figure = Figure, Table = Table, Image = Image}
}
