library(extrafont)
# font_import() ## need to run on new system to make database of fonts. 
# fonts() ## just useful to see all of the fonts that I have. 
# font_import(pattern ="trebuc", prompt=FALSE) # can't find this font anymore. Come back to this later. 
# loadfonts()
theme_JAG_presentation <- function (base_size = 12){
  theme_grey(base_size=base_size) %+replace%
    theme(
      rect = element_rect(colour = "black", fill = NA, size = 0.5, linetype = 1, inherit.blank = TRUE),
      text = element_text(family = "Helvetica",
                          face = "plain",
                          color = "black",
                          size = base_size,
                          hjust = 0.5,
                          vjust = 0.5, 
                          angle = 0, 
                          lineheight = 0.9,
                          margin = margin(), 
                          debug = FALSE),
      line = element_line(colour = "black",
                        size = 0.5,
                        linetype = 1,
                        lineend = "butt"),
      #axis.line = element_blank(), 
      axis.line = element_line(colour = 'grey', size = 0.3),
      axis.ticks = element_line(colour = "grey", size = 0.2), 
      axis.title.x = element_text(size = base_size+2, colour = 'black', vjust = 0), 
      axis.title.y = element_text(face="bold", size=base_size+2, colour = 'black', angle = 90, vjust = 0.5), 
      axis.text.x = element_text(size = base_size, colour = "black"),
      axis.text.y = element_text(size = base_size, colour = "black"),
      legend.background = element_rect(colour = NA, fill = 'white'), 
      legend.key = element_rect(colour = NA, fill = 'white'), 
      legend.key.height = NULL, 
      legend.key.width = NULL,     
      legend.text = element_text(colour = 'black'), 
      legend.title = element_text(face = "bold", hjust = 0, colour = 'black'), 
      legend.position = "right", 
      legend.text.align = NULL, 
      legend.title.align = NULL, 
      legend.direction = "vertical", 
      legend.box = NULL,    
      
     # panel.background = element_rect(fill = "white") ,
     panel.background = element_rect(fill = NA) , ## to keep points from being hidden
  #    panel.background = element_rect(fill = "black", colour = "white"), 
  #    panel.border = element_rect(fill = NA, colour = "white"), 
      panel.border = element_rect(fill = NA, colour = "white"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      
      
      strip.background = element_rect(fill = "grey30", colour = "grey10"), 
      strip.text.x = element_text(colour = 'black'), 
      strip.text.y = element_text(colour = 'black', angle = -90), 
      
     ## plot.background = element_rect(colour = 'white', fill = 'white'), 
  plot.background = element_rect(colour = NA, fill = NA), # to avoid covering points
      plot.title = element_text(colour = "black")
    )
}
