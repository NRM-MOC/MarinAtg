library(tidyverse)
library(knitr)
library(kableExtra)
library(shinydashboard)
source("functions.R")

data <- read_csv("aggdata.csv", guess_max = 50000)

ui <- dashboardPage(
  header = dashboardHeader(),
  # Input
  sidebar = dashboardSidebar(
    selectInput("Art", "Art", choices = data %>% pull(art) %>% unique()),
    selectInput("Substans", "Variabel", choices = data %>% arrange(var_order) %>% pull(label) %>% unique()),
    sliderInput("Start", "Start", min = 2000, max = 2012, value = 2005, sep = "")
  ),
  body <-  dashboardBody(
    fluidRow(
      tabBox(width = 12, title =  htmlOutput("title"),
             tabPanel("Table", downloadButton("downloadTable", "Download"), htmlOutput("table")),
             tabPanel("Plot", plotOutput("plot", width = "800px", height = "600px")),
             tabPanel("Power plot", plotOutput("powerplot", width = "600px", height = "600px"))
      ))))

server <- function(input, output){
  download_table <- reactive({
    table_data() %>%
      get_table() %>% 
      select(basin,
             Location = LOC, 
             `Yearly change (%)` = relative_rate, 
             `Confidence interval (95%)` = CI,
             `Predicted value 2018` = predicted_2018,
             `p-value equal slopes` = p_vs_sat, 
             `Power vs 10% inc.` = power,
             `Limit` = limit,
             `Time to limit (years)` = ttt,
             `n obs` = n,
             `n all < LOD` = n.all.lod)
  })
  table_data <- reactive({
    data %>% 
      filter(label == input[["Substans"]], 
             art == input[["Art"]],
             YEAR >= input[["Start"]])
  })
  power_table <- reactive({
    tab <- data %>% 
      filter(label == input[["Substans"]], 
             art == input[["Art"]],
             YEAR >= 2000) %>%     
      group_by(basin) %>% 
      get_table() %>% 
      select(LOC, basin, power, var) %>% 
      mutate(YEAR = 2000)
    for (year in 2001:2012) {
      tab <- bind_rows(tab, data %>% 
                         filter(label == input[["Substans"]], 
                                art == input[["Art"]],
                                YEAR >= year) %>%     
                         group_by(basin) %>% 
                         get_table() %>% 
                         select(LOC, basin, power, var) %>% 
                         mutate(YEAR = year)
      )
    }
    tab
  }) 
  output[["table"]] <- renderText(
    ifelse(nrow(table_data()) == 0, "Not available", 
           table_data() %>%     
             get_table() %>% 
             arrange(basin_order) %>% 
             select(basin,
                    Location = LOC, 
                    `Yearly change (%)` = relative_rate, 
                    `Confidence interval (95%)` = CI,
                    `Predicted value 2018` = predicted_2018,
                    `p-value equal slopes` = p_vs_sat, 
                    `Power vs 10% inc.` = power,
                    `Limit` = limit,
                    `Time to limit (years)` = ttt,
                    `n obs` = n,
                    `n all < LOD` = n.all.lod) %>% 
             pretty_basin_kable(caption = NULL))
  )
  output[["plot"]] <- renderPlot(
    try(
      table_data() %>%     
        get_table() %>% 
        mutate(func = pmap(list(intercept, slope, first_year), ~tibble(YEAR = ..3:2018, conc = exp(..1 + ..2*(..3:2018))))) %>% 
        unnest(func) %>% 
        ggplot(aes(x = YEAR, y = conc, group = LOC)) +       
        geom_point(data = table_data(), aes(x = YEAR, y = value), size = .1, color = "darkgrey") +
        geom_line(aes(linetype = (basin != LOC)), color = "steelblue") +
        theme_minimal() + facet_wrap(~basin) + 
        theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5),
              panel.grid.minor = element_blank()) + 
        scale_x_continuous(breaks = input[["Start"]]:2018, limits = c(input[["Start"]], 2018)) +
        xlab("") + ggtitle("Fitted concentration curves (Basins (solid), Stations (dashed))"),
      TRUE)
  )
  output[["powerplot"]] <- renderPlot(
    try(
      power_table() %>% filter(basin == LOC) %>% 
        ggplot(aes(x = YEAR, y = power, color = LOC))  + theme_minimal() + 
        geom_point() + geom_line() + 
        theme(legend.position = "top", legend.title = element_blank(),
              panel.grid.minor = element_blank()) + 
        scale_x_continuous(breaks = 2000:2012) + 
        xlab("Starting year") + ylab("Power against 10% increase") + 
        ylim(c(0, 1)),
      TRUE)
  )
  output[["title"]] <- renderText(
    paste0(input[["Art"]], " (", mocis_get_unit_HTML(filter(data, label == input[["Substans"]]) %>% pull(var) %>% .[1], input[["Art"]]), ")")
  )
  output[["downloadTable"]] <- downloadHandler(
    filename = paste0(input[["Substans"]],"_in_", input[["Art"]], "_from_", input[["Start"]], ".csv"),
    content = function(file) {
      write.csv2(download_table(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)