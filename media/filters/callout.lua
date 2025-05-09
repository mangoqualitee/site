function Callout(el)
  if quarto.doc.isFormat("html") then
    -- Set collapse to true only if the callout has content and collapse isn't already set
    -- if el.content and not el.collapse then
    if not el.collapse then
      el.collapse = true
    end
    return el
  end
end