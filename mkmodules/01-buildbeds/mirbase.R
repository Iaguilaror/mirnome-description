# cargamos pacman
library( "pacman" )

# cargamos libs
p_load( "vroom",
        "dplyr",
        "tidyr",
        "stringr" )

## Read args from command line
args = commandArgs( trailingOnly = TRUE )

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "*.gff3" ## the input gff3

## Passing args to named objects
input_gff <- args[1]

# colnames
nombres_columnas <- c( "chromosome",
                       "source",
                       "feature",
                       "start",
                       "end",
                       "score",
                       "strand",
                       "frame",
                       "attribute" )

# cargamos datos
whole <- vroom( file = input_gff,
                delim = "\t",
                comment = "#",
                col_names = nombres_columnas )

# separate the annotation
separated <- separate( data = whole,
                       col = attribute,
                       into = c( "ID",
                                 "Alias",
                                 "Name",
                                 "Derives" ),
                       sep = ";" ) %>%
  select( -source, -frame ) %>%
  mutate( ID = str_remove( string = ID,
                           pattern = "ID=" ),
          Alias = str_remove( string = Alias,
                              pattern = "Alias=" ),
          Name = str_remove( string = Name,
                             pattern = "Name=" ),
          Derives = str_remove( string = Derives,
                                pattern = "Derives_from=" ) ) %>%
  mutate( Derives = ifelse( test = is.na(Derives),
                            yes = ID,
                            no = Derives ) )

# operate on upstream and downstram regiones
primaries <- separated %>%
  filter( feature == "miRNA_primary_transcript" )

# define plus
surrounding_plus <- primaries %>%
  filter( strand == "+" ) %>%
  mutate ( up1k_start = start - 1e3,
           up1k_end = start - 1 ) %>%
  mutate( down1k_start = end + 1,
          down1k_end = end + 1e3 )

# define down
surrounding_down <- primaries %>%
  filter( strand == "-" ) %>%
  mutate ( up1k_start = end + 1,
           up1k_end = end + 1e3 ) %>%
  mutate( down1k_start = start - 1e3,
          down1k_end = start - 1 )

# gather
surrounding <- rbind( surrounding_down,
                      surrounding_plus ) %>%
  mutate( order = case_when( chromosome == "chr1" ~ 1,
                             chromosome == "chr2" ~ 2,
                             chromosome == "chr3" ~ 3,
                             chromosome == "chr4" ~ 4,
                             chromosome == "chr5" ~ 5,
                             chromosome == "chr6" ~ 6,
                             chromosome == "chr7" ~ 7,
                             chromosome == "chr8" ~ 8,
                             chromosome == "chr9" ~ 9,
                             chromosome == "chr10" ~ 10,
                             chromosome == "chr11" ~ 11,
                             chromosome == "chr12" ~ 12,
                             chromosome == "chr13" ~ 13,
                             chromosome == "chr14" ~ 14,
                             chromosome == "chr15" ~ 15,
                             chromosome == "chr16" ~ 16,
                             chromosome == "chr17" ~ 17,
                             chromosome == "chr18" ~ 18,
                             chromosome == "chr19" ~ 19,
                             chromosome == "chr20" ~ 20,
                             chromosome == "chr21" ~ 21,
                             chromosome == "chr22" ~ 22,
                             chromosome == "chrX" ~ 23,
                             chromosome == "chrY" ~ 24,
                             chromosome == "chrMT" ~ 25) )

# build a bed for up1k
bed_up <- surrounding %>%
  select( chromosome,
          up1k_start,
          up1k_end,
          Name,
          score,
          strand,
          order ) %>%
  mutate( Name = paste0( Name, "_upstream_1k") ) %>%
  arrange( order, up1k_start, up1k_end  ) %>%
  select( -order )

# build a bed for down1k
bed_down <- surrounding %>%
  select( chromosome,
          down1k_start,
          down1k_end,
          Name,
          score,
          strand,
          order ) %>%
  mutate( Name = paste0( Name, "_downstream_1k") ) %>%
  arrange( order, down1k_start, down1k_end  ) %>%
  select( -order )

# save beds
write.table( x = bed_up,
             file = "up1k_primary_mirBase_22.bed",
             quote = FALSE, sep = "\t", row.names = FALSE,
             col.names = FALSE )

# save beds
write.table( x = bed_down,
             file = "down1k_primary_mirBase_22.bed",
             quote = FALSE, sep = "\t", row.names = FALSE,
             col.names = FALSE )

#### get primary coordinates only ----
bed_primaries <- separated %>%
  filter( feature == "miRNA_primary_transcript" ) %>%
  select( chromosome,
          start,
          end,
          Name,
          score,
          strand ) %>%
  mutate( Name = paste0( Name, "_primary" ) )

# save beds
write.table( x = bed_primaries,
             file = "primary_mirBase_22.bed",
             quote = FALSE, sep = "\t", row.names = FALSE,
             col.names = FALSE )

### create a bed for every position of a mature mirna ----
mature_tmp <- separated %>%
  filter( feature == "miRNA" )

# first plus
mature_plus_tmp <- mature_tmp %>%
  filter( strand == "+" )

# create empty df
mature_untangled_tmp_plus <- mature_plus_tmp[0,] %>%
  mutate( mirposition = as.numeric(  ) )

# start loop
for ( row_in_turn in 1:nrow( mature_plus_tmp) ) {

  message( "untangling ", mature_plus_tmp$ID[ row_in_turn ] )

  # sacar el rango
  position.v <- mature_plus_tmp$start[ row_in_turn ]:mature_plus_tmp$end[ row_in_turn ]

  # longitud
  longitud <- length( position.v )

  # sacar el cromosoma
  chrom.v <- rep( mature_plus_tmp$chromosome[ row_in_turn ], longitud )

  # sacar el feature
  feature.v <- rep( mature_plus_tmp$feature[ row_in_turn ], longitud )

  # sacar el score
  score.v <- rep( ".", longitud )

  # sacar el strand
  strand.v <- rep( mature_plus_tmp$strand[ row_in_turn ], longitud )

  # sacar el ID
  id.v <- rep( mature_plus_tmp$ID[ row_in_turn ], longitud )

  # sacar el Alias
  alias.v <- rep( mature_plus_tmp$Alias[ row_in_turn ], longitud )

  # sacar el Name
  name.v <- rep( mature_plus_tmp$Name[ row_in_turn ], longitud )

  # sacar el Derives
  derives.v <- rep( mature_plus_tmp$Derives[ row_in_turn ], longitud )

  # create an intermediate df
  df_tmp <- data.frame( chromosome = chrom.v,
                        feature = feature.v,
                        start = position.v -1,
                        end = position.v,
                        score = score.v,
                        strand = strand.v,
                        ID = id.v,
                        Alias = alias.v,
                        Name = name.v,
                        Derives = derives.v,
                        mirposition = 1:longitud )

  # bind intermediate df
  mature_untangled_tmp_plus <- rbind( mature_untangled_tmp_plus,
                                      df_tmp )
}

### Now with minus strand -----
# now minus
mature_minus_tmp <- mature_tmp %>%
  filter( strand == "-" )

# create empty df
mature_untangled_tmp_minus <- mature_minus_tmp[0,] %>%
  mutate( mirposition = as.numeric(  ) )

# start loop
for ( row_in_turn in 1:nrow( mature_minus_tmp ) ) {

  message( "untangling ", mature_minus_tmp$ID[ row_in_turn ] )

  # sacar el rango
  position.v <- mature_minus_tmp$start[ row_in_turn ]:mature_minus_tmp$end[ row_in_turn ]

  # longitud
  longitud <- length( position.v )

  # sacar el cromosoma
  chrom.v <- rep( mature_minus_tmp$chromosome[ row_in_turn ], longitud )

  # sacar el feature
  feature.v <- rep( mature_minus_tmp$feature[ row_in_turn ], longitud )

  # sacar el score
  score.v <- rep( ".", longitud )

  # sacar el strand
  strand.v <- rep( mature_minus_tmp$strand[ row_in_turn ], longitud )

  # sacar el ID
  id.v <- rep( mature_minus_tmp$ID[ row_in_turn ], longitud )

  # sacar el Alias
  alias.v <- rep( mature_minus_tmp$Alias[ row_in_turn ], longitud )

  # sacar el Name
  name.v <- rep( mature_minus_tmp$Name[ row_in_turn ], longitud )

  # sacar el Derives
  derives.v <- rep( mature_minus_tmp$Derives[ row_in_turn ], longitud )

  # create an intermediate df
  df_tmp <- data.frame( chromosome = chrom.v,
                        feature = feature.v,
                        start = position.v -1,
                        end = position.v,
                        score = score.v,
                        strand = strand.v,
                        ID = id.v,
                        Alias = alias.v,
                        Name = name.v,
                        Derives = derives.v,
                        mirposition = rev( 1:longitud ) )

  # bind intermediate df
  mature_untangled_tmp_minus <- rbind( mature_untangled_tmp_minus,
                                       df_tmp )
}

## gather plus and minus ----
allmature <- bind_rows( mature_untangled_tmp_plus,
                        mature_untangled_tmp_minus ) %>%
  mutate( order = case_when( chromosome == "chr1" ~ 1,
                             chromosome == "chr2" ~ 2,
                             chromosome == "chr3" ~ 3,
                             chromosome == "chr4" ~ 4,
                             chromosome == "chr5" ~ 5,
                             chromosome == "chr6" ~ 6,
                             chromosome == "chr7" ~ 7,
                             chromosome == "chr8" ~ 8,
                             chromosome == "chr9" ~ 9,
                             chromosome == "chr10" ~ 10,
                             chromosome == "chr11" ~ 11,
                             chromosome == "chr12" ~ 12,
                             chromosome == "chr13" ~ 13,
                             chromosome == "chr14" ~ 14,
                             chromosome == "chr15" ~ 15,
                             chromosome == "chr16" ~ 16,
                             chromosome == "chr17" ~ 17,
                             chromosome == "chr18" ~ 18,
                             chromosome == "chr19" ~ 19,
                             chromosome == "chr20" ~ 20,
                             chromosome == "chr21" ~ 21,
                             chromosome == "chr22" ~ 22,
                             chromosome == "chrX" ~ 23,
                             chromosome == "chrY" ~ 24,
                             chromosome == "chrMT" ~ 25) ) %>%
  arrange( order, start, end  ) %>%
  select( -order ) %>%
  mutate( region = ifelse( test = mirposition <= 8 &
                             mirposition >= 2 ,
                           yes = "seed",
                           no = "mature" ) ) %>%
  mutate( Name = paste( Name,
                        region,
                        mirposition,
                        sep = "_" ) )

# prepare bedfile
#### get primary coordinates only ----
bed_matures <- allmature %>%
  select( chromosome,
          start,
          end,
          Name,
          score,
          strand )

# save beds
write.table( x = bed_matures,
             file = "matures_perbase_mirBase_22.bed",
             quote = FALSE, sep = "\t", row.names = FALSE,
             col.names = FALSE )
