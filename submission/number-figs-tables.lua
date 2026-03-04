-- Lua filter for pandoc: adds Figure N / Table N numbering to captions
-- Tracks counters for figures and tables sequentially

local fig_counter = 0
local tab_counter = 0

function Figure(el)
  fig_counter = fig_counter + 1
  -- Prepend "Figure N. " to the caption
  if el.caption and el.caption.long and #el.caption.long > 0 then
    local first_block = el.caption.long[1]
    if first_block and first_block.content then
      table.insert(first_block.content, 1, pandoc.Strong(pandoc.Str("Figure " .. fig_counter .. ". ")))
    end
  end
  return el
end

function Table(el)
  tab_counter = tab_counter + 1
  -- Prepend "Table N. " to the caption
  if el.caption and el.caption.long and #el.caption.long > 0 then
    local first_block = el.caption.long[1]
    if first_block and first_block.content then
      table.insert(first_block.content, 1, pandoc.Strong(pandoc.Str("Table " .. tab_counter .. ". ")))
    end
  end
  return el
end
