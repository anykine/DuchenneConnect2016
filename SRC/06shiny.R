# shiny app

library(shiny)

shinyUI(
  fluidPage(
    titlePanel("title"),
    
    sidebarLayout(
      sidebarPanel("side panel"),
      mainPanel("main")
    )
  )
)