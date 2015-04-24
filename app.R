#!/usr/bin/env Rscript

#' #Introduction
#' Purpose of this is to identity KOs with disproportiate expression given its gDNA background

library(MetamapsDB)
library(tidyr)
library(ggplot2)
library(ggvis)
library(shiny)
library(dplyr)
theme_set(theme_bw())

drawGradient = function(values, colors) {
    mvalue = median(values, na.rm=T)

    lowerMid = colorRampPalette(colors[1:2])
    midUpper= colorRampPalette(colors[2:3])
    #20 colors from btm to top; spread be taken

    lowerV = cut(values[which(values < mvalue)], breaks=10)
    levels(lowerV) = lowerMid(10)

    upperV = cut(values[which(values > mvalue)], breaks=10)
    levels(upperV) =  midUpper(10)

    colorVector = values
    colorVector[which(values < mvalue)] = as.character(lowerV)
    colorVector[which(values > mvalue)] = as.character(upperV)


    colorVector[is.na(colorVector)] = "#000000"
    colorVector
}

#' ## Top500 KOs
#' Loads list of top500 kos & convert into ggvis obj

#data(top500kos)
#edgeDF = ig2ggvis(prettifyGraph(grepgraph(top500kos))) #only if you're on water
load("edgeDF.rda")

ggplot_edge = edgeDF %>% select(x:y, row)
ge = setNames(
              do.call(rbind,lapply(unique(ggplot_edge$row), function(rowID){
                                   df = subset(ggplot_edge, row == rowID)
                                   cbind(df[1,1:2], df[2,1:2])
            })),
              c("x", "y","xend", "yend"))

pointsDF = edgeDF %>% select(x,y, name, label, type) %>% unique
koDF = edgeDF     %>% select(x,y, name, label, type) %>% unique %>% filter(type == 'ko')
cpdDF = edgeDF    %>% select(x,y, name, label, type) %>% unique %>% filter(type=='cpd')
#' ## gDNA and mRNA information
#load("out/graph_prototype.rda")
load("kodf.rda")

#all = do.call(rbind,lapply(data, function(eachKO) eachKO$df))
#all %>% filter(readType == 'mRNA')

app = shinyApp(
####################################################################################################
    ui = fluidPage(
        #Title #################################################
        titlePanel("Top500KOs"),

        #Input #################################################
        sidebarLayout(
            sidebarPanel(
                        sliderInput('perc', h4('Select NXX:'), value = 0.9, min = 0.1, max = 0.95, step=0.05),
                        sliderInput('ratio', 'gDNA to cDNA Ratio:', value=0, min=0, max=3, step=0.05),
                        sliderInput("range", "Show KO labels within ratio range:", min = 0, max = 4, value = c(0.8,1.8), step=0.05), width=2
                    ),
        #Plots #################################################
        mainPanel(
              tabsetPanel(
                  tabPanel("gDNA X cDNA",
                           p("Ratio"),
                           fluidRow(
                                    column(12,plotOutput("NXXGVM")),
                                    column(12,plotOutput("Ratio"))
                                    )
                           ),
                tabPanel("Top500",
                         fluidRow(
                    #column(12,plotOutput("graph", width=800, height=500)),
                    column(12, ggvisOutput("top500"))
                         )
                         )
                  )
              )
        )
        ),
####################################################################################################
        server = function(input, output){
            all_values <- function(x) {
                if (is.null(x)) 
                    return(NULL)
                paste0(names(x), ": ", format(x), collapse = "<br />")
            }

            Ndf = reactive({
                tt = all                                       %>%
                #all %>%
                    filter(percentage == input$perc)      %>%
                    select(contigsRequired, ko, readType) %>%
                    spread(readType, contigsRequired)     %>%
                    mutate(ratio = cDNA/gDNA ,
                            ko   = gsub("^", "ko:", ko))  %>%
                    filter(ratio > input$ratio)
            })

#Plot1: Raw num scatter
            output$NXXGVM = renderPlot({
                ggplot(Ndf(), aes(gDNA, cDNA)) +
                geom_point()
            })
#Plot2: Ratio
            output$Ratio = renderPlot({
                ggplot(Ndf(), aes(x=as.factor(1), y=ratio))                                                     +
                    geom_jitter(aes(color=ratio))                                                               +
                    xlab("KOs")                                                                                 +
                    ggtitle(sprintf("Ratio: cDNA / gDNA (%s)", input$perc))                                     +
                    scale_color_gradient2("Presence", midpoint=1, low="#91bfdb", mid="#ffffbf", high="#fc8d59")
            })

#Plot3: GGVIS
            reactive({
                #' keeps KOs in the metabolic graph
                koDF = merge(koDF, #kos in the graph
                             #tt %>% select(ko, ratio),
                             Ndf() %>% select(ko, ratio),
                             by.x="name", by.y="ko", all.x=T)
                koDF$included = 1; koDF$included[is.na(koDF$ratio)] = 0 #may not have the contig data for these KOs
                #koDF$fill = drawGradient(koDF$ratio, c("#e41a1c", "#377eb8", "#4daf4a"))
                #koDF$fill = drawGradient(koDF$ratio, rev(c("#fc8d59", "#ffffbf", "#91bfdb")))
                koDF$fill = drawGradient(koDF$ratio, rev(c("#d73027", "#ffffbf", "#4575b4"))) #red blue
                
                #cat(koDF$ratio)
                edgeDF %>%
                    ggvis(~x, ~y)                                                                                                                                                               %>%
                    group_by(row)                                                                                                                                                           %>%
                    layer_paths()                                                                                                                                                           %>%
                    #ko
                    #layer_points(data = koDF, size:= 50, shape= ~factor(included)) %>%
                    layer_points(data = koDF, size:= 50, shape= ~factor(included), fill:=~fill, stroke:="black") %>% #, fill= ~ratio)                                                                                                                                %>%
                    #cpd
                    layer_points(data = cpdDF, size:= 10, fill := "grey")                                                                                                                                %>%
                    layer_text(data=koDF, text:=~ratio, fontSize:=1) %>%
                    layer_text(data=koDF %>% filter(ratio > input$range[1] & ratio < input$range[2]),text:=~label, fontSize:=10) %>%
                    add_axis("x",title = "", properties = axis_props(axis = list(strokeWidth = 0),grid = list(strokeWidth = 0), ticks = list(strokeWidth = 0),labels = list(fontSize = 0))) %>%
                    add_axis("y",title = "", properties = axis_props(axis = list(strokeWidth = 0),grid = list(strokeWidth = 0), ticks = list(strokeWidth = 0),labels = list(fontSize = 0))) %>%
                    #scale_ordinal("fill", range = c("blue","red"))                                                                                                                          %>%
                    set_options(height = 800, width = 800) %>%
                    add_tooltip(all_values, "hover")
            }) %>% bind_shiny("top500")

#Plot4: Graph
#            output$graph= renderPlot({
#                koDF = merge(koDF, Ndf() %>% select(ko, ratio), by.x="name", by.y="ko", all.x=T)
#                
#                #when the KO is not included in the analysis (but cant be, just in case)
#                koDF$included = 1; koDF$included[is.na(koDF$ratio)] = 0
#
#                ggplot()                                                                                        +
#                    geom_segment(data = ge, aes(x=x, xend=xend, y=y, yend=yend))                                +
#                    geom_point(data = koDF, size=6, aes(x=x, y=y, color=ratio, shape=as.factor(included)))      +
#                    geom_point(data = cpdDF, size= 2, fill = "grey", aes(x=x, y=y))                             +
#                    scale_color_gradient2("Presence", midpoint=1, low="#91bfdb", mid="#ffffbf", high="#fc8d59") +
#                    theme(line = element_blank(),
#                          title = element_blank())
#            })
        })
