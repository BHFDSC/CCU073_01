######################################################################################
## TITLE: Create supplementary figure 1 -Flow diagram for cohrt selection process
##
## Author: Genevieve Cezard
######################################################################################

install.packages(c("DiagrammeRsvg", "rsvg"))
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

setwd("2. Projects/2. CVD-COVID-UK/CCU073 project/Code")

g <- grViz("
digraph population_flow {
  graph [layout = dot, rankdir = TB]

  node [shape = box, style = filled, fillcolor = white, fontname = Helvetica]

  A [label = 'People with a record in GDPPR on January 1st 2022\n(n = 65,494,390)']
  B [label = 'With sex, date of birth and LSOA*\n and living in England\n(n = 65,414,515)']
  C [label = 'Passed quality assurance**\n(n = 65,396,840)']
  D [label = 'People aged 40 years old and over\n(n = 30,628,405)']
  E [label = 'People with no CVD event prior baseline\n(n = 24,733,400)']

  A -> B
  B -> C
  C -> D
  D -> E
}
")
# Convert to SVG then PNG
svg <- export_svg(g)
rsvg_png(charToRaw(svg), file = "CCU073_Flow_diagram_v3.png")