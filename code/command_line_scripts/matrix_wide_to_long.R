suppressMessages(require(data.table))
suppressMessages(require(stringr))
suppressMessages(require(argparser))

#print(.libPaths())


parser <- arg_parser('Convert a square matrix into a long one', hide.opts = TRUE)

parser <- add_argument(parser, "--cores", type="integer", default=NULL,
                       help="Indicate how many cores should be used for computation. If not set, data.table reads environment variables and uses all ligcal CPUs available")
parser <- add_argument(parser, "--input", nargs=Inf,
                       help='Specify the path to at least one square matrix or use wildcard \"*\" to process multiple chromosomes')
parser <- add_argument(parser, "--outfile", default=getwd(),
                       help="For processing of multiple files: Indicate path to which output tables should be saved \n 
                       If only a single file is processed, you can sepcify an path + filename")
parser <- add_argument(parser, "--distanceCutoff",type="integer", 
                       help='Specify a maximum distance of considered interactions in bp. Must be at least twice the value of the maximally considered IS distance! E.g. To remove contacts exceeding 2Mb of distance, specify -d 2000000')
parser <- add_argument(parser, "--expectedContacts", default='cis',
                       help='If the supplied matrix is genome wide (and contains interchromosomal contacts), please indicate "trans".')


main <- function() {
  args <- parse_args(parser)
  filenames <- args$input
  scope <- args$expectedContacts
  
  #catch error: scope is neither cis nor trans
  if(sum(scope != c('cis', 'trans')) > 1){
    stop('expectedContacts argument must be either "cis" or "trans"')
  }  
  
  #catch error: no matrix supplied
  if(is.na(filenames[1])){
    stop('NO input matrices supplied to function. Please specify at least one matrix in long format for calculation of IS ')
  }
  
  #Print out how many matrices will be converted
  if(length(filenames)==1){
    print('1 matrix will be converted to long format')
  }
  if(length(filenames) > 1){
    print(paste(length(filenames), 'matrices will be converted to long format'))
  }
  
  ###Convert cis matrices
  if(scope == 'cis'){
    distanceCutoff = args$distanceCutoff
    
    #Print if distance cutoff was supplied. Only use cutoff for cis-matrices
    if(! is.na(distanceCutoff)){
      print(paste('Removing contacts exceeding', distanceCutoff, 'basepairs'))
    }
    
    
    #loop through all specified files
    for (filename in filenames) {
      print(paste('read in file:',filename))
      #have to do it the unpractical way as often the matrices contain multiple NAs and the colClasses are estimated as 'character' otherwise
      colnames_wide <- colnames(fread(file = filename, nrows = 1))
      mat_wide <- fread(file = filename, colClasses = list(double=2:length(colnames_wide)))
      
      #catch error: input matrix not 'square'
      if(dim(mat_wide)[1]+1 != dim(mat_wide)[2]){
        stop('Input matrix does not meet shape requirements. Please input a matrix that has coordinates (chr:start-end) as first column so that it has 1 more column than rows ')
      }
      
      #reshape matrix and split location identifiers in 3 separate columns
      mat_long <-melt(mat_wide, id.vars = colnames_wide[1], measure = colnames_wide[-1])
      mat_long[, c('A_chrom', 'A_start','A_end') := tstrsplit(get(colnames_wide[1]), ':|-', type.convert=TRUE)]
      mat_long[, c('B_chrom', 'B_start','B_end') := tstrsplit(variable, ':|-', type.convert=TRUE)]
      mat_long[, distance := abs(A_start-B_start)]
  
      #catch error: trans chromosomal matrix
      if(length(unique(mat_long$A_chrom)) != 1 |
         length(unique(mat_long$B_chrom)) != 1){
        stop('matrix has trans chromosomal interactions. Reshaping to long format with 3 column output (start, stop, value) would loose information')
      }
      
      #Determine which chromosome
      chromosome_ID_A <- unique(mat_long$A_chrom)
      chromosome_ID_B <- unique(mat_long$B_chrom)
      #catch error: more than 1 chromosome in matrix
      if(length(chromosome_ID_A) > 1|
         length(chromosome_ID_B) > 1){
        stop('matrix column or rownames contain more than one chromosome identifer')
      }
      #catch error: no cis interactions
      if(chromosome_ID_A != chromosome_ID_B){
        stop('Chromosome identiiers of rows and columns not matching')
      }
      
      #remove bins above distance threshold (if specified with args$distanceCutoff)
      if(! is.na(distanceCutoff)){
        mat_long <- mat_long[distance <= distanceCutoff]
      }
      
      
      #select A_start, B_start, value
      mat_long_three_column <- mat_long[,c('A_start','B_start', 'value')]
      #mat_long_three_column <- mat_long_three_column[!is.na(mat_long_three_column$value)]#drop empty values
      mat_long_three_column <- mat_long_three_column[mat_long_three_column$A_start != mat_long_three_column$B_start] #remove diagonal
      
      
      #generate filename:
      filename_without_path <- unlist(lapply(str_split(filename, pattern='/'), tail, 1))
      #filename_without_fileending <- unlist(lapply(str_split(filename_without_path, pattern='\\.'), head, 1))
      filename_long_format <- paste(filename_without_path, 'long.tsv.gz', sep="_")
      if(length(filenames)==1){
        filename_outpath = args$outfile
      }
      if(length(filenames)>1){
        filename_outpath = paste0(args$outfile, filename_long_format)
      }
      
      #write to output directory
      #catch error: no writing permission in output directory
      if(file.create(filename_outpath) != TRUE){
        stop('no writing permission for the output directory')
      }
      fwrite(mat_long_three_column, filename_outpath, sep = '\t', na=NA)
    }
  }
  
  ###Convert trans matrices
    if(scope == 'trans'){
      
      #generate filename:
      filename_outpath = args$outfile
      #catch error: no writing permission in output directory
      if(file.create(filename_outpath) != TRUE){
        stop('no writing permission for the output directory')
      }
      
      #loop through all specified files
      for (filename in filenames) {
      print(paste('read in file:',filename))
      colnames_wide <- colnames(fread(file = filename, nrows = 1))
      mat_wide <- fread(file = filename, colClasses = list(double=2:length(colnames_wide)))
      
      #catch error: input matrix not 'square'
      if(dim(mat_wide)[1]+1 != dim(mat_wide)[2]){
        stop('Input matrix does not meet shape requirements. Please input a matrix that has coordinates (chr:start-end) as first column so that it has 1 more column than rows ')
      }
      
      #reshape matrix and split location identifiers in 3 separate columns
      mat_long <- melt(mat_wide, id.vars = colnames_wide[1], measure = colnames_wide[-1])
      mat_long[, c('A_chrom', 'A_start','A_end') := tstrsplit(get(colnames_wide[1]), ':|-', type.convert=TRUE)]
      mat_long[, c('B_chrom', 'B_start','B_end') := tstrsplit(variable, ':|-', type.convert=TRUE)]
      mat_long[,'interaction_type' := fcase(A_chrom == B_chrom, 'cis',
                                            default = 'trans')]

      #select A_chrom, A_start, B_chrom, B_start, interaction_type, value
      mat_long_three_column <- mat_long[,c('A_chrom', 'A_start', 'B_chrom', 'B_start', 'interaction_type', 'value')]
      #remove NAs
      mat_long_three_column <- mat_long_three_column[!is.na(mat_long_three_column$value)]
      #remove diagonal 
      mat_long_cis <- mat_long_three_column[mat_long_three_column$A_chrom == mat_long_three_column$B_chrom]
      mat_long_cis <- mat_long_cis[mat_long_cis$A_start != mat_long_cis$B_start] 
      #gl
      mat_long_three_column <- rbindlist(list(mat_long_cis, 
                                         mat_long_three_column[mat_long_three_column$A_chrom != mat_long_three_column$B_chrom]))
      #sort table
      mat_long_three_column <- mat_long_three_column[order(A_chrom, A_start, B_chrom, B_start)]
      
      #write to output directory
      fwrite(mat_long_three_column, filename_outpath, sep = '\t', na=NA)
    }
  }
}

main()
