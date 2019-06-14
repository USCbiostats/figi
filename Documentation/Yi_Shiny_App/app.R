#setwd('H:/DataHarm/Applications/ShinyR/FIGI')
#setwd("/home/ylin2/DataHarm/Applications/ShinyR/FIGI/")
rm(list = ls())

#setwd("/Users/mak/Dropbox/n_figi/FIGI Epi data harmonization shiny app update/")
setwd("/Users/mak/Dropbox/FIGI/Documentation/Yi_Shiny_App/")
load('Initial_figi.Rdata')
library(shiny)
library(ggplot2)
library(rmeta)
library(xtable)

ui <- fluidPage(
  tabsetPanel(
    tabPanel('Counts',
      headerPanel('FIGI Epidemiologic Data Exploration (June 2017)'),
      sidebarPanel(
        uiOutput("Box1"),
        uiOutput("Box2")),
      mainPanel(
        uiOutput('counts')
      )),
    tabPanel('Boxplot-p1',plotOutput("boxplot1",height = "650px")),
    tabPanel('Boxplot-p2',plotOutput("boxplot2",height = "650px")),
    tabPanel('forestplot',plotOutput("forestplot",height = "750px")),
    tabPanel('Funnelplot',plotOutput("funnelplot",height = "650px"))
  )
)

server <- function(input, output) {

  output$Box1 = renderUI(selectInput('Categ1', 'Category',unique(var.ls$Group),'Study'))
  output$Box2 = renderUI({
    if(is.null(input$Categ1)){
      return()
    }else{
       selectInput('Categ2', 'Variable',
             var.ls$Description[var.ls$Group %in% input$Categ1],'Case/Control status')}
    })

  output$counts <- renderUI({
    if(is.null(input$Categ2) | is.null(input$Categ1)){return()
      }else{ if(input$Categ2 %in% ''){return()
            }else{
                  if(var.ls$Short.Name[var.ls$Description %in% input$Categ2] %in% c(NA,NULL,var.sub)){
                   plot(0,0,type="n",axes=F,xlab='',ylab='',main='')
                  }else{
                        M<-Counts[[var.ls$Short.Name[var.ls$Description %in% input$Categ2]]]

                        if(var.ls$Short.Name[var.ls$Description %in% input$Categ2] %in% 'outc'){
                            M <- print(xtable(M,align=c(rep("l|",3),rep("r", ncol(M)-2))),
                                      include.rownames=F,floating=FALSE, tabular.environment="array",
                                      comment=FALSE, print.results=FALSE, hline.after=c(-1,0,nrow(M)-1,nrow(M)))#,
                        }else{
                            M <- print(xtable(M,align=c(rep("l|",3),rep("c", (ncol(M)-1)/2-1),
                                        "c|",rep("c", (ncol(M)-1)/2))),include.rownames=F,include.colnames=F,
			                                  floating=FALSE, tabular.environment="array",comment=FALSE,
			                                 print.results=FALSE,hline.after=c(-1,0,1,2,nrow(M)))
                        }
                        html <- paste0("$$", M, "$$")
                        withMathJax(HTML(html))
                }
            }}})


  output$boxplot1 <- renderPlot({
    if(is.null(input$Categ2) | is.null(input$Categ1)){return()
    }else{ if(input$Categ2 %in% ''){return()
    }else{
      if(var.ls$Short.Name[var.ls$Description %in% input$Categ2] %in% c(NA,NULL,var.sub)){
        plot(0,0,type="n",axes=F,xlab='',ylab='',main='')
        text(0,0,'Not available for adjusted variables\nOr look at results for * subgroup',
             cex=4)
      }else{
            Box1[[var.ls$Short.Name[var.ls$Description %in% input$Categ2]]]
    }}}})

  output$boxplot2 <- renderPlot({
    if(is.null(input$Categ2) | is.null(input$Categ1)){return()
    }else{ if(input$Categ2 %in% ''){return()
    }else{
      if(var.ls$Short.Name[var.ls$Description %in% input$Categ2] %in% c(NA,NULL,var.sub)){
        plot(0,0,type="n",axes=F,xlab='',ylab='',main='')
        text(0,0,'Not available for adjusted variables\nOr look at results for * subgroup',
             cex=4)
      }else{
        Box2[[var.ls$Short.Name[var.ls$Description %in% input$Categ2]]]
  }}}})

  output$forestplot <- renderPlot({
    if(is.null(input$Categ2) | is.null(input$Categ1)){return()
    }else{ if(input$Categ2 %in% ''){return()
    }else{
      if(var.ls$Short.Name[var.ls$Description %in% input$Categ2] %in%
         setdiff(var.ls$Short.Name,meta$var)){
        plot(0,0,type="n",axes=F,xlab='',ylab='',main='')
        text(0,0,'Not available for adjusted variables\nOr look at results for * subgroup',
             cex=4)
      }else{
          plot.forestplot(meta$var[meta$Description %in% input$Categ2],
                          meta.std=meta,study=res.tot,lab=lab,lab.std=lab.std,sample=sample)
      }
  }}})

  output$funnelplot <- renderPlot({
    if(is.null(input$Categ2) | is.null(input$Categ1)){return()
    }else{ if(input$Categ2 %in% ''){return()
    }else{
      if(var.ls$Short.Name[var.ls$Description %in% input$Categ2] %in%
         setdiff(var.ls$Short.Name,meta$var)){
        plot(0,0,type="n",axes=F,xlab='',ylab='',main='')
        text(0,0,'Not available for adjusted variables\nOr look at results for * subgroup',
             cex=4)
      }else{
      plot.funnel(meta$var[meta$Description %in% input$Categ2],dat=data.fun)
  }}}})

}

shinyApp(ui = ui, server = server)
