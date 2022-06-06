# cargamos pacman
library( "pacman" )

# cargamos libs
p_load( "vroom",
        "dplyr",
        "tidyr",
        "stringr",
        "scales",
        "cowplot",
        "ggplot2" )

## Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "test/data" ## the input dir

## Passing args to named objects
bed_dir <- args[1]

# create a function to read each file
allfiles <- list.files( path = bed_dir,
                        pattern = ".\\.bed",
                        full.names = TRUE )

basedf <- vroom( allfiles,
                 col_names = FALSE ) %>% 
  rename( "chr" = 1,
          "start" = 2,
          "end" = 3,
          "feature" = 4,
          "score" = 5,
          "strand" = 6,
          "chr2" = 7,
          "start2" = 8,
          "end2" = 9,
          "variant" = 10,
          "score2" = 11,
          "strand2" = 12,
          "overlap" = 13 )

#
filtered <- basedf %>% 
  filter( overlap != 0 ) %>% 
  select( -score,
          -chr2,
          -start2,
          -end2,
          -score2,
          -strand2,
          -overlap )

# separated
separated <- separate( data = filtered,
                       col = feature,
                       into = c( "mir", "feature", "featposition" ),
                       sep = "\\." ) %>% 
  separate( data = .,
            col = variant,
            into = c( "var", "AF" ),
            sep = "_AF=" ) %>% 
  separate( data = .,
            col = var,
            into = c( "varchr", "varpos", "ref", "alt" ),
            sep = "_",
            remove = FALSE ) %>% 
  mutate( featposition = as.numeric( featposition ),
          varpos = as.numeric( varpos ),
          AF = as.numeric( AF) )

# count number of regs
varcount <- separated %>% 
  group_by( var ) %>% 
  summarise( count = n( ) )

# uniqvars
uniqvars <- varcount %>% 
  filter( count == 1 ) %>% 
  pull( var )

# multivars
multivars <- varcount %>% 
  filter( count > 1 ) %>% 
  pull( var )

# find uniq regs
uniqdf <- separated %>% 
  filter( var %in% uniqvars )

# find dup regs
dupdf <- separated %>% 
  filter( var %in% multivars )


#remove upstream or downstream because it means a  primary or seed is present
filtereddupdf <- dupdf %>% 
  mutate( importance = case_when( feature == "seed" ~ 1,
                                  feature == "mature" ~ 2,
                                  feature == "primary" ~ 3,
                                  feature == "upstream_1k" ~ 4,
                                  feature == "downstream_1k" ~ 5 ) ) %>% 
  arrange( var, importance ) %>% 
  group_by( var ) %>%
  slice( 1 ) %>% 
  select( -importance ) %>% 
  ungroup( )

# gather uniqs and deduped
cleaned <- bind_rows( uniqdf,
                      filtereddupdf ) %>% 
  mutate( feature = factor( x = feature, levels = c( "upstream_1k",
                                                     "primary",
                                                     "mature",
                                                     "seed",
                                                     "downstream_1k" ) ) )

# recount to make sure no dups are found
varcount2 <- cleaned %>% 
  group_by( var ) %>% 
  summarise( count = n( ) )

# save the final file
write.table( x = cleaned,
             file = "allvariants.bed",
             quote = FALSE, sep = "\t",
             row.names = FALSE, col.names = TRUE)

### preview frequencies
box1 <- ggplot( data = cleaned,
                mapping = aes( x = feature,
                               y = AF ) ) +
  geom_boxplot( outlier.size = 0.5,
                # size = 0.5,
                width = 0.5,
                color = "gray30"
  ) +
  scale_x_discrete( limits = c( "upstream_1k",
                                "primary",
                                "mature",
                                "seed",
                                "downstream_1k" ) ) +
  scale_y_continuous( limits = c(0, 1),
                      breaks = seq( from = 0, to = 1, by =0.25 ),
                      labels = percent ) +
  labs( title = "miRNA Population Allele Frequencies",
        x = "Sequence type",
        y = "Allele Frequency" ) +
  theme_classic( )

# try a jitter
tagged <- cleaned %>% 
  mutate( interest = ifelse( test = (feature == "mature" | feature == "seed") & AF > 0.05,
                             yes = "voi",
                             no = "non_voi") )

# try a jittered plot
boxjitter1 <- ggplot( data = filter( tagged, interest == "non_voi"),
                      mapping = aes( x = feature,
                                     y = AF ) ) +
  geom_hline( yintercept = 0.05,
              linetype = "dashed",
              color = "tomato" ) +
  geom_point( size = 1, shape = 20,
              alpha = 0.4,
              position = position_jitter( seed = 77 ) ) +
  geom_point( data = filter(tagged, interest == "voi"),
              size = 2, shape = 21,
              alpha = 0.7, fill = "tomato",
              position = position_jitter( seed = 77 ) ) +
  scale_x_discrete( limits = c( "upstream_1k",
                                "primary",
                                "mature",
                                "seed",
                                "downstream_1k" ) ) +
  scale_y_continuous( limits = c(0, 1.01),
                      breaks = seq( from = 0, to = 1, by =0.25 ),
                      labels = percent ) +
  labs( title = "miRNA Population Allele Frequencies",
        x = "Sequence type",
        y = "Allele Frequency" ) +
  theme_bw( ) +
  theme( axis.text.x = element_blank( ),
         axis.ticks.x = element_blank( ),
         panel.grid.major.x = element_blank(),
         legend.position = "none" ) +
  facet_wrap( ~ feature, nrow = 1, strip.position = "bottom" )

# try a panel
panel1 <- plot_grid( box1, boxjitter1, ncol = 1 )

ggsave( plot = panel1,
        filename = "boxplots.svg",
        width = 10,
        height = 7 )

# summarise vois
tagged %>% 
  group_by( feature, interest ) %>% 
  summarise( count = n( ) ) %>% 
  write.csv( x = ., quote = FALSE,
             file = "variants_of_interest_by_type.csv",
             row.names = FALSE )
