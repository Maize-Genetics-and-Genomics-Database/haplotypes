library(shiny)
library(feather)


hmp_samples = feather("./hapblock_info_summary_new_hier.feather") 
hmp_samples = unname(unlist(unique(hmp_samples[,"sample"])))


fluidPage(
tags$head(includeCSS("/srv/shiny-server/dropdown.css")),
    includeCSS("./styles.css"),
    tabsetPanel(
        tabPanel("Single Sample View",
                 headerPanel("Haplotype Assignments by Sample"),
                 sidebarLayout(
                     sidebarPanel(
                         fluidRow(
                             selectInput('hmp_sample', 'Sample', sort(hmp_samples), width = "175px", selectize = F, size = 15),
                             br(),
                             selectInput('scale_select', 'Scale', c("Genetic","Physical"), width = "175px")
                         ),
                         width = 2
                         
                     ),
                     mainPanel(
                         fluidRow(   
                             conditionalPanel(condition = "input.scale_select == 'Genetic'", plotOutput('plot')),
                             conditionalPanel(condition = "input.scale_select == 'Physical'", plotOutput('plot_phys'))
                         )
                     ),
                     position = "left")
        ),
        tabPanel("Multiple Samples View",
                 headerPanel("Comparison of Haplotype Assignments Across Multiple Samples"),
                 sidebarLayout(
                     sidebarPanel(
                         fluidRow(
                             selectizeInput('comp_inbreds', 'Select Samples to Compare', sort(hmp_samples), multiple = T, width = "175px"),
                             br(),
                             selectInput('chrom_select', 'Chromosome', seq(1,10,1), width = "175px"),
                             br(),
                             selectInput('scale_select2', 'Scale', c("Genetic","Physical"), width = "175px"),
                             br(),
                             checkboxInput("checkbox2", "Show All Chroms"),
                             br(),
                             actionButton('show_multi_plot',"Plot", width = "175px")
                             
                         ),
                         width = 2
                         
                     ),
                     mainPanel(
                         fluidRow(   
                             conditionalPanel(condition = "input.checkbox2 == false & input.scale_select2 == 'Genetic'", uiOutput('plot2.ui')),
                             conditionalPanel(condition = "input.checkbox2 == false & input.scale_select2 == 'Physical'", uiOutput('plot2_phys.ui')),
                             conditionalPanel(condition = "input.checkbox2 == true & input.scale_select2 == 'Genetic'", plotOutput('plot3')),
                             conditionalPanel(condition = "input.checkbox2 == true & input.scale_select2 == 'Physical'", plotOutput('plot3_phys'))
                         )
                     ),
                     position = "left")
        ),
        tabPanel("Compare Samples to Reference",
                 headerPanel("Comparison of Haplotype Assignments to a Reference Sample"),
                 sidebarLayout(
                     sidebarPanel(
                         fluidRow(
                             selectInput('ref_inbred', 'Select Reference Sample', sort(hmp_samples), width = "175px"),
                             selectizeInput('comp_inbreds2', 'Samples to Compare', sort(hmp_samples), multiple = T, width = "175px"),
                             br(),
                             selectInput('chrom_select2', 'Chromosome', seq(1,10,1), width = "175px"),
                             br(),
                             checkboxInput("checkbox3", "Show All Chroms"),
                             br(),
                             selectInput('scale_select3', 'Scale', c("Genetic","Physical"), width = "175px"),
                             br(),
                             actionButton('show_multi_plot2',"Plot", width = "175px")
                         ),
                         width = 2
                         
                     ),
                     mainPanel(
                         fluidRow(   
                             conditionalPanel(condition = "input.checkbox3 == false & input.scale_select3 == 'Genetic'", uiOutput('plot4.ui')),
                             conditionalPanel(condition = "input.checkbox3 == false & input.scale_select3 == 'Physical'", uiOutput('plot4_phys.ui')),
                             conditionalPanel(condition = "input.checkbox3 == true & input.scale_select3 == 'Genetic'", plotOutput('plot5')),
                             conditionalPanel(condition = "input.checkbox3 == true & input.scale_select3 == 'Physical'", plotOutput('plot5_phys'))
                         )
                     ),
                     position = "left")
        ))
    

)
