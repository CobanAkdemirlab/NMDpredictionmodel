library(ggplot2)
library(dplyr)

cm <- cm %>%
    mutate(
        # White text on dark/high-percentage cells,
        # black text on light/low-percentage cells
        text_color = ifelse(Percentage >= 50, "white", "black")
    )

ggplot(cm, aes(x = Predicted, y = True, fill = Percentage)) +
    
    geom_tile(
        color = "white",
        linewidth = 2
    ) +
    
    geom_text(
        aes(
            label = Label,
            color = text_color
        ),
        size = 4,
        fontface = "bold",
        lineheight = 1.2
    ) +
    
    scale_color_identity() +
    
    # LOW = light, HIGH = dark
    scale_fill_gradient(
        name = "Proportion (%)",
        low  = "#DCEAF7",
        high = "#17365D",
        limits = c(25, 75),
        breaks = c(30, 40, 50, 60, 70),
        guide = guide_colorbar(
            reverse = FALSE
        )
    ) +
    
    labs(
        title = "Confusion Matrix @ 0.42",
        x = "Predicted",
        y = "True"
    ) +
    
    coord_fixed() +
    
    theme_classic(base_size = 14) +
    
    theme(
        plot.title = element_text(
            hjust = 0.5,
            face = "bold"
        ),
        axis.title = element_text(
            face = "bold"
        ),
        axis.text = element_text(
            size = 14,
            face = "bold"
        ),
        legend.title = element_text(
            face = "bold",
            size = 10
        )
    )
