library(shiny)
library(feather)
library(ggplot2)
library(data.table)
library(dplyr)
library(scales)

plottheme = theme(axis.text.y = element_text(size=12, color="#292929"), axis.text.x = element_text(size=12, color="#292929"),
                  axis.title.x = element_text(size=14, face="bold"), axis.title.y = element_text(size=14, face="bold"),
                  plot.title = element_text(size=14, face="bold"), strip.text.x = element_text(size=10, face="bold"), 
                  strip.background = element_rect(colour="#A2B5CD", fill="white"),
                  panel.background = element_rect(colour="black", fill="white"),
                  legend.text = element_text(size = 12), legend.background = element_rect(fill=alpha('white', 0.5), colour = 'black'),
                  legend.title = element_text(size = 12), legend.justification=c(0, 1.7),
                  panel.border = element_rect(colour = "#D3D3D3", fill = NA, size = 0.6),
                  plot.background = element_rect(colour = "black", size = 1))

function(input, output) {

## IMPORT DATA
    
    hmp_inbred_data = read_feather("./hapblock_info_summary_new_hier.feather")
    sample_metadata = fread("./sample_metadata.csv")
    hapblock_pos = fread("./map_positions.csv", header = T)
    hmp_inbred_data = left_join(hmp_inbred_data, sample_metadata, by = "sample")
    hmp_inbred_data = left_join(hmp_inbred_data, hapblock_pos[,c("chr","startpos","genpos_start")], by = c("chr","startpos"))
    hmp_inbred_data = left_join(hmp_inbred_data, hapblock_pos[,c("chr","endpos","genpos_end")], by = c("chr","endpos"))
    maxgenpos = max(hmp_inbred_data$genpos_end)
    maxphyspos = max(hmp_inbred_data$endpos)/1000000
    
    print(head(hmp_inbred_data))
     
## PREP INBRED DATA
    inbred_data <- reactive({
        # inbred = "282set_Mt42"
        inbred = input$hmp_sample
        a = hmp_inbred_data[which(hmp_inbred_data$sample == inbred),]
        a = a[with(a, order(chr,block)),]
        print(head(a))
        a$Hex_Code = as.character(a$Hex_Code)
        a$origin_sample = as.character(a$origin_sample)
        # founderhap_levels = data.table(table(a[,c("origin_sample","Hex_Code")]))
        founderhap_levels = a %>%
            group_by(origin_sample, Hex_Code) %>%
            summarise(N = sum((endpos/1000000)-(startpos/1000000)))
        founderhap_levels = if(length(which(founderhap_levels$N == 0)) > 0){
            founderhap_levels[-which(founderhap_levels$N == 0),]
        } else {founderhap_levels}
        founderhap_levels = founderhap_levels[with(founderhap_levels, order(-N)),]
        a$origin_sample = factor(a$origin_sample, levels = founderhap_levels$origin_sample)
        a$Hex_Code = factor(a$Hex_Code, levels = founderhap_levels$Hex_Code)
        return(a)
    })
    
    multi_inbred_data <- reactive({
        input$show_multi_plot
        
        isolate({
        if(length(input$comp_inbreds) > 0){
            
            # inbreds = c("LH74","B73","A632")
            inbreds = input$comp_inbreds
            print(inbreds)
            a = hmp_inbred_data[which(hmp_inbred_data$sample %in% inbreds),]
            print(head(a))
            a = a[with(a, order(chr,block)),]
            a$Hex_Code = as.character(a$Hex_Code)
            a$origin_sample = as.character(a$origin_sample)
            founderhap_levels = a %>%
                group_by(origin_sample, Hex_Code) %>%
                summarise(N = sum((endpos/1000000)-(startpos/1000000)))
            founderhap_levels = if(length(which(founderhap_levels$N == 0)) > 0){
                founderhap_levels[-which(founderhap_levels$N == 0),]
            } else {founderhap_levels}
            founderhap_levels = founderhap_levels[with(founderhap_levels, order(-N)),]
            a$origin_sample = factor(a$origin_sample, levels = founderhap_levels$origin_sample)
            a$Hex_Code = factor(a$Hex_Code, levels = founderhap_levels$Hex_Code)
            # a$sample = factor(a$sample, levels = sort(inbreds)) # alphabetic sort
            a$sample = factor(a$sample, levels = inbreds) # sort based on input sample order
            print(head(a))
            print(str(a))
            return(a)
        } else {
            NULL
        }
        })
    })
    
    
    multi_inbred_ref_data <- reactive({
        if(length(input$ref_inbred) > 0 & length(input$comp_inbreds2) > 0){
            # inbreds = c("PHG39","207","A632")
            # ref_inbred = "B14"
            inbreds = input$comp_inbreds2
            
            ref <- hmp_inbred_data[which(hmp_inbred_data$sample == input$ref_inbred),]
            # ref <- hmp_inbred_data[which(hmp_inbred_data$sample == "B14"),]
            ref = if(length(which(duplicated(ref$startpos))) > 0) {
                ref[-which(duplicated(ref$startpos)),]
            } else {ref}
            
            ref_hex_code = unique(ref$Hex_Code)
            if(length(ref_hex_code) != 1) {print("Something went wrong when getting the ref_inbred hex code")} else {}
            
            b <- hmp_inbred_data[which(hmp_inbred_data$sample %in% inbreds),]
            # b <- hmp_inbred_data[which(hmp_inbred_data$sample %in% c("PHG39","PH207","A632")),]
            
            b = rbind(ref,b)
            
            # the idea here is that if the comparison inbreds have the same origin inbred as the reference inbred, then they match to B73
            # so anywhere where they are equal, the origin sample for the comparison sample at that chrom and block can be changed to the name of the ref inbred
            # the hex code also needs to change
            colnames(ref)[which(colnames(ref) == "origin_sample")] <- "ref_origin_sample"
            colnames(ref)[which(colnames(ref) == "sample")] <- "ref_inbred"
            
            a = full_join(b, ref[,c("chr","block","ref_inbred","ref_origin_sample")], by = c("chr", "block"))
            
            same_haps = which(a$ref_origin_sample == a$origin_sample)
            diff_haps = which((a$ref_origin_sample != a$origin_sample) | (is.na(a$ref_origin_sample) & !(is.na(a$origin_sample))))
            
            # modify values in place
            # a[same_haps,"origin_sample"] <- "B14"
            a[same_haps,"origin_sample"] <- input$ref_inbred
            a$Hex_Code = as.character(a$Hex_Code)
            a[same_haps,"Hex_Code"] <- "#00CD00"
            
            diff_hex = "#FBEC5D"
            a[diff_haps,"Hex_Code"] <- diff_hex

            df_dups <- a[,c("chr","startpos", "sample")]
            
            a = if(length(which(duplicated(df_dups))) > 0) {
                a[-which(duplicated(df_dups)),]
            } else {a}
            
            # # add back in the reference sample at the top so it shows up first in the plot; make colors == ref sample
            # # ref$ref_origin_sample <- "B73"
            # ref$ref_origin_sample <- input$ref_sample
            # ref$Hex_Code <- ref_hex_code
            # colnames(ref)[which(colnames(ref) == "ref_origin_sample")] <- "origin_sample"
            # colnames(ref)[which(colnames(ref) == "ref_inbred")] <- "sample"
            # a = rbind.fill(ref,a)
            
            a = a[with(a, order(chr,block)),]
            a$Hex_Code = as.character(a$Hex_Code)
            a$origin_sample = as.character(a$origin_sample)
            founderhap_levels = data.table(table(a[,c("origin_sample","Hex_Code")]))
            founderhap_levels = if(length(which(founderhap_levels$N == 0)) > 0){
                founderhap_levels[-which(founderhap_levels$N == 0),]
            } else {founderhap_levels}
            founderhap_levels = founderhap_levels[with(founderhap_levels, order(-N)),]
            a$origin_sample = factor(a$origin_sample, levels = founderhap_levels$origin_sample)
            
            a$Hex_Code = factor(a$Hex_Code, levels = c("#00CD00","#FBEC5D"))
            # a$sample = factor(a$sample, levels = c("B14","PHG39","PH207","A632"))
            # a$sample = factor(a$sample, levels = c(input$ref_inbred,sort(input$comp_inbreds2))) # alphabetic sort
            a$sample = factor(a$sample, levels = c(input$ref_inbred, inbreds)) # sort based on input sample order
            print("Multi Inbred Ref Data")
            print(head(a))
            # print(str(a))
            return(a)
        } else {
            NULL
        }
    })
    

    
    
## MAKE PLOTS
    
## PLOT TO VIEW ONE INBRED, ALL CHROMS - GENETIC
    output$plot <- renderPlot({

        p <- ggplot(inbred_data()) +
            geom_rect(aes(xmin = genpos_start, xmax = genpos_end, fill = Hex_Code, ymin = 0, ymax = 1)) +
            scale_fill_identity(guide = "legend", labels = levels(inbred_data()$origin_sample)[1:15],
                                breaks = levels(inbred_data()$Hex_Code)[1:15]) +
            coord_flip() +
            scale_x_reverse(limits = c(maxgenpos+10, 0), breaks = seq(0,maxgenpos+10, 10)) +
            xlab("Genetic Position (cM)") +
            ylab("") +
            facet_wrap(~chr, nrow = 1) +
            plottheme +
            ggtitle(paste("Inbred:", input$hmp_sample, "  Pedigree:", unique(inbred_data()$pedigree), "   ", unique(inbred_data()$comment))) +
            theme(legend.position = c(0, 0),
                  plot.margin = unit(c(25, 25, 105, 25), "points"),
                  axis.text.x = element_blank(), axis.ticks = element_blank()) +
            guides(fill = guide_legend(title = "Hapgroups", nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1))
        
        
        print(p)

        
        }, height = 700, width = 1100)
    

## PLOT TO VIEW ONE INBRED, ALL CHROMS - PHYSICAL
    output$plot_phys <- renderPlot({

        p <- ggplot(inbred_data()) +
            geom_rect(aes(xmin = startpos/1000000, xmax = endpos/1000000, fill = Hex_Code, ymin = 0, ymax = 1)) +
            scale_fill_identity(guide = "legend", labels = levels(inbred_data()$origin_sample)[1:15],
                                breaks = levels(inbred_data()$Hex_Code)[1:15]) +
            coord_flip() +
            scale_x_reverse(limits = c(maxphyspos+10,0), breaks = seq(0,maxphyspos+50, 50)) +
            xlab("Physical Position (Mb)") +
            ylab("") +
            facet_wrap(~chr, nrow = 1) +
            plottheme +
            ggtitle(paste("Inbred:", input$hmp_sample, "  Pedigree:", unique(inbred_data()$pedigree), "   ", unique(inbred_data()$comment))) +
            theme(legend.position = c(0, 0),
                  plot.margin = unit(c(25, 25, 105, 25), "points"),
                  axis.text.x = element_blank(), axis.ticks = element_blank()) +
            guides(fill = guide_legend(title = "Hapgroups", nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1))
        
        
        print(p)

        
            }, height = 700, width = 1100)
    
    
    
    
    comp_inb_plot_width <- reactive({
        input$show_multi_plot
        
        isolate({
        a = 200+(length(unique(multi_inbred_data()$sample))*100)
        return(a)
        })
    })
    
    
## PLOT TO VIEW MULTIPLE INBREDS, SINGLE CHROM, GENETIC
    
    output$plot2 <- renderPlot({
        
        input$show_multi_plot
        
        isolate({
            if(length(input$comp_inbreds) > 0 & input$checkbox2 == FALSE & input$scale_select2 == "Genetic"){

        dd <- multi_inbred_data()
        dd <- dd[which(dd$chr == input$chrom_select),]
                
        
        p <- ggplot(dd) +
            geom_rect(aes(xmin = genpos_start, xmax = genpos_end, fill = Hex_Code, ymin = 0, ymax = 1)) +
            scale_fill_identity(guide = "legend", labels = levels(dd$origin_sample)[1:15],
                                breaks = levels(dd$Hex_Code)[1:15]) +
            coord_flip() +
            scale_x_reverse(limits = c(maxgenpos+10, 0), breaks = seq(0,maxgenpos+10, 10)) +
            xlab("Genetic Position (cM)") +
            ylab("") +
            facet_wrap(~sample, nrow = 1) +
            plottheme +
            theme(legend.position = c(1.05, 1.2),
                  plot.margin = unit(c(25, 155, 25, 25), "points"),
                  axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(size = 16, angle = 90)) +
            guides(fill = guide_legend(title = "Hapgroups", ncol = 1, label.position = "right", 
                                       keywidth = 2, keyheight = 1, direction = "vertical"))
        
        
        print(p)
        
        } else {
            NULL
        }
    })
    })
    
    
    output$plot2.ui <- renderUI({
        input$show_multi_plot
        
        isolate({
            
        plotOutput("plot2", height = 800, width = comp_inb_plot_width())

        })
        
    })
    

    ## PLOT TO VIEW MULTIPLE INBREDS, SINGLE CHROM, PHYSICAL
    
    output$plot2_phys <- renderPlot({
        
        input$show_multi_plot
        
        isolate({
            if(length(input$comp_inbreds) > 0 & input$checkbox2 == FALSE & input$scale_select2 == "Physical"){
                
                dd <- multi_inbred_data()
                dd <- dd[which(dd$chr == input$chrom_select),]
                
                p <- ggplot(dd) +
                    geom_rect(aes(xmin = startpos/1000000, xmax = endpos/1000000, fill = Hex_Code, ymin = 0, ymax = 1)) +
                    scale_fill_identity(guide = "legend", labels = levels(dd$origin_sample)[1:15],
                                        breaks = levels(dd$Hex_Code)[1:15]) +
                    coord_flip() +
                    scale_x_reverse(limits = c(maxphyspos+10,0), breaks = seq(0,maxphyspos+50, 50)) +
                    xlab("Physical Position (Mb)") +
                    ylab("") +
                    facet_wrap(~sample, nrow = 1) +
                    plottheme +
                    theme(legend.position = c(1.05, 1.2),
                          plot.margin = unit(c(25, 155, 25, 25), "points"),
                          axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(size = 16, angle = 90)) +
                    guides(fill = guide_legend(title = "Hapgroups", ncol = 1, label.position = "right", 
                                               keywidth = 2, keyheight = 1, direction = "vertical"))
                
                
                print(p)
                
            } else {
                NULL
            }
        })
    })
    
    
    output$plot2_phys.ui <- renderUI({
        input$show_multi_plot
        
        isolate({
            
            plotOutput("plot2_phys", height = 800, width = comp_inb_plot_width())

        })
        
    })
    
    
    
## PLOT TO VIEW MULTIPLE INBREDS, ALL CHROMS, GENETIC
    output$plot3 <- renderPlot({
        input$show_multi_plot
        
        isolate({
        if(length(input$comp_inbreds) > 0 & input$checkbox2 == TRUE & input$scale_select2 == "Genetic"){
            
            p <- ggplot(data = multi_inbred_data()) +
                geom_rect(aes(xmin = genpos_start, xmax = genpos_end, fill = Hex_Code, ymin = as.integer(sample), ymax = as.integer(sample)+0.9)) +
                scale_fill_identity(guide = "legend", labels = levels(multi_inbred_data()$origin_sample)[1:15],
                                    breaks = levels(multi_inbred_data()$Hex_Code)[1:15]) +
                coord_flip() +
                scale_x_reverse(limits = c(maxgenpos+10, 0), breaks = seq(0,maxgenpos+10, 10)) +
                scale_y_continuous(breaks = 1:length(levels(multi_inbred_data()$sample)), labels = levels(multi_inbred_data()$sample), 
                                   position = "right") +
                xlab("Genetic Position (cM)") +
                ylab("") +
                facet_wrap(~chr, nrow = 2, strip.position = "left") +
                plottheme +
                theme(legend.position = c(1.05, 1.2),
                      plot.margin = unit(c(25, 175, 25, 25), "points"),
                      axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(size = 16, angle = 90)) +
                guides(fill = guide_legend(title = "Hapgroups", ncol = 1, label.position = "right", 
                                           keywidth = 2, keyheight = 1, direction = "vertical"))
            
                
                print(p)
        } else {
            NULL
        }
    })}, height = 800, width = 1300)
    
    ## PLOT TO VIEW MULTIPLE INBREDS, ALL CHROMS, PHYSICAL
    output$plot3_phys <- renderPlot({
        input$show_multi_plot
        
        isolate({
            if(length(input$comp_inbreds) > 0 & input$checkbox2 == TRUE & input$scale_select2 == "Physical"){
                
                p <- ggplot(data = multi_inbred_data()) +
                    geom_rect(aes(xmin = startpos/1000000, xmax = endpos/1000000, fill = Hex_Code, ymin = as.integer(sample), ymax = as.integer(sample)+0.9)) +
                    scale_fill_identity(guide = "legend", labels = levels(multi_inbred_data()$origin_sample)[1:15],
                                        breaks = levels(multi_inbred_data()$Hex_Code)[1:15]) +
                    coord_flip() +
                    scale_x_reverse(limits = c(maxphyspos+10,0), breaks = seq(0,maxphyspos+50, 50)) +
                    scale_y_continuous(breaks = 1:length(levels(multi_inbred_data()$sample)), labels = levels(multi_inbred_data()$sample), 
                                       position = "right") +
                    xlab("Physical Position (Mb)") +
                    ylab("") +
                    facet_wrap(~chr, nrow = 2, strip.position = "left") +
                    plottheme +
                    theme(legend.position = c(1.05, 1.2),
                          plot.margin = unit(c(25, 175, 25, 25), "points"),
                          axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(size = 16, angle = 90)) +
                    guides(fill = guide_legend(title = "Hapgroups", ncol = 1, label.position = "right", 
                                               keywidth = 2, keyheight = 1, direction = "vertical"))
                
                
                print(p)
            } else {
                NULL
            }
        })}, height = 800, width = 1300)
    
    

    
## PLOT TO VIEW MULTIPLE INBREDS USING A REFERENCE INBRED, SINGLE CHROM, GENETIC
    comp_inb2_plot_width <- reactive({
        input$show_multi_plot2
        
        isolate({
            
        a = 200+(length(unique(input$comp_inbreds2))*100)
        return(a)
        
        })
    })
    
    output$plot4 <- renderPlot({
        input$show_multi_plot2
        
        isolate({
        if(length(input$ref_inbred) > 0 & length(input$comp_inbreds2) > 0 & input$checkbox3 == FALSE & input$scale_select3 == "Genetic"){

            dd <- multi_inbred_ref_data()
            dd <- dd[which(dd$chr == input$chrom_select2),]
            
            p <- ggplot(dd) +
                geom_rect(aes(xmin = genpos_start, xmax = genpos_end, fill = Hex_Code, ymin = 0, ymax = 1)) +
                scale_fill_identity(guide = "legend", labels = c("Same Haplotype","Different Haplotype"),
                                    breaks = levels(multi_inbred_ref_data()$Hex_Code)[1:2]) +
                coord_flip() +
                scale_x_reverse(limits = c(maxgenpos+10, 0), breaks = seq(0,maxgenpos+10, 10)) +
                xlab("Genetic Position (cM)") +
                ylab("") +
                facet_wrap(~sample, nrow = 1) +
                plottheme +
                ggtitle(paste("Reference Inbred:", input$ref_inbred)) +
                theme(legend.position = c(0, 0),
                      plot.margin = unit(c(25, 25, 105, 25), "points"),
                      axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(size = 16, angle = 90)) +
                guides(fill = guide_legend(title = "", nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1))
            
            
            print(p)
        
        } else {
            NULL
        }
        })
        
    })
    
    
    output$plot4.ui <- renderUI({
        input$show_multi_plot2
        
        isolate({
        if(length(input$ref_inbred) > 0 & length(input$comp_inbreds2) > 0 & input$checkbox3 == FALSE & input$scale_select3 == "Genetic"){
        
        plotOutput("plot4", height = 800, width = comp_inb2_plot_width())
        
        } else { NULL }
        })
    })   
    
    
    
## PLOT TO VIEW MULTIPLE INBREDS USING A REFERENCE INBRED, SINGLE CHROM, PHYSICAL
    output$plot4_phys <- renderPlot({
        input$show_multi_plot2
        
        isolate({
            if(length(input$ref_inbred) > 0 & length(input$comp_inbreds2) > 0 & input$checkbox3 == FALSE & input$scale_select3 == "Physical"){
                
                
                dd <- multi_inbred_ref_data()
                dd <- dd[which(dd$chr == input$chrom_select2),]
                
                p <- ggplot(dd) +
                    geom_rect(aes(xmin = startpos/1000000, xmax = endpos/1000000, fill = Hex_Code, ymin = 0, ymax = 1)) +
                    scale_fill_identity(guide = "legend", labels = c("Same Haplotype","Different Haplotype"),
                                        breaks = levels(multi_inbred_ref_data()$Hex_Code)[1:2]) +
                    coord_flip() +
                    scale_x_reverse(limits = c(maxphyspos+10,0), breaks = seq(0,maxphyspos+50, 50)) +
                    xlab("Physical Position (Mb)") +
                    ylab("") +
                    facet_wrap(~sample, nrow = 1) +
                    plottheme +
                    ggtitle(paste("Reference Inbred:", input$ref_inbred)) +
                    theme(legend.position = c(0, 0),
                          plot.margin = unit(c(25, 25, 105, 25), "points"),
                          axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(size = 16, angle = 90)) +
                    guides(fill = guide_legend(title = "", nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1))
                
                
                print(p)
                
            } else {
                NULL
            }
        })
        
    })
    
    
    output$plot4_phys.ui <- renderUI({
        input$show_multi_plot2
        
        isolate({
            if(length(input$ref_inbred) > 0 & length(input$comp_inbreds2) > 0 & input$checkbox3 == FALSE & input$scale_select3 == "Physical"){
                
                plotOutput("plot4_phys", height = 800, width = comp_inb2_plot_width())
                
            } else { NULL }
        })
    })   
    
    
    
## PLOT TO VIEW MULTIPLE INBREDS USING A REFERENCE INBRED, ALL CHROMS, GENETIC
    output$plot5 <- renderPlot({
        input$show_multi_plot2
        
        isolate({
        if(length(input$ref_inbred) > 0 & length(input$comp_inbreds2) > 0 & input$checkbox3 == TRUE & input$scale_select3 == "Genetic"){

               p <- ggplot(data = multi_inbred_ref_data()) +
                    geom_rect(aes(xmin = genpos_start, xmax = genpos_end, fill = Hex_Code, ymin = as.integer(sample), ymax = as.integer(sample)+0.9)) +
                    scale_fill_identity(guide = "legend", labels = c("Same Haplotype","Different Haplotype"),
                                       breaks = levels(multi_inbred_ref_data()$Hex_Code)[1:2]) +
                    coord_flip() +
                    scale_x_reverse(limits = c(maxgenpos+10, 0), breaks = seq(0,maxgenpos+10, 10)) +
                    scale_y_continuous(breaks = 1:length(levels(multi_inbred_ref_data()$sample)), labels = levels(multi_inbred_ref_data()$sample), 
                                       position = "right") +
                    xlab("Genetic Position (cM)") +
                    ylab("") +
                    facet_wrap(~chr, nrow = 2, strip.position = "left") +
                    plottheme +
                    ggtitle(paste("Reference Inbred:", input$ref_inbred)) +
                    theme(legend.position = c(0, 0),
                          plot.margin = unit(c(25, 25, 105, 25), "points"),
                          axis.text.x.top = element_text(size = 14, angle = 90, vjust = 1, hjust = 0),
                          axis.ticks=element_blank(), strip.text.x = element_text(size = 16),
                          strip.text.y = element_text(size = 12, angle = -180)) +
                    guides(fill = guide_legend(title = "", nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1))

            print(p)
            
            } else {
                NULL
            }
        })
        }, height = 800, width = 1300)

    ## PLOT TO VIEW MULTIPLE INBREDS USING A REFERENCE INBRED, ALL CHROMS, PHYSICAL
    output$plot5_phys <- renderPlot({
        input$show_multi_plot2
        
        isolate({
            if(length(input$ref_inbred) > 0 & length(input$comp_inbreds2) > 0 & input$checkbox3 == TRUE & input$scale_select3 == "Physical"){
                
                p <- ggplot(data = multi_inbred_ref_data()) +
                    geom_rect(aes(xmin = startpos/1000000, xmax = endpos/1000000, fill = Hex_Code, ymin = as.integer(sample), ymax = as.integer(sample)+0.9)) +
                    scale_fill_identity(guide = "legend", labels = c("Same Haplotype","Different Haplotype"),
                                        breaks = levels(multi_inbred_ref_data()$Hex_Code)[1:2]) +
                    coord_flip() +
                    scale_x_reverse(limits = c(maxphyspos+10,0), breaks = seq(0,maxphyspos+50, 50)) +
                    scale_y_continuous(breaks = 1:length(levels(multi_inbred_ref_data()$sample)), labels = levels(multi_inbred_ref_data()$sample), 
                                       position = "right") +
                    xlab("Physical Position (Mb)") +
                    ylab("") +
                    facet_wrap(~chr, nrow = 2, strip.position = "left") +
                    plottheme +
                    ggtitle(paste("Reference Inbred:", input$ref_inbred)) +
                    theme(legend.position = c(0, 0),
                          plot.margin = unit(c(25, 25, 105, 25), "points"),
                          axis.text.x.top = element_text(size = 14, angle = 90, vjust = 1, hjust = 0),
                          axis.ticks=element_blank(), strip.text.x = element_text(size = 16),
                          strip.text.y = element_text(size = 12, angle = -180)) +
                    guides(fill = guide_legend(title = "", nrow = 1, label.position = "bottom", keywidth = 2, keyheight = 1))
                
                print(p)
                
            } else {
                NULL
            }
        })
    }, height = 800, width = 1300)
    
}
