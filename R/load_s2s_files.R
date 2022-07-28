#' generate output dir
#'
#' This function generates two dirs, one 'analysis_output' and within this a directory
#' with the current date
#' @param output_dir_location Path to the location you want the directories to be generated
#' @return the path to the directory with the current data
#' @export
generate_output_dir <- function(output_dir_location) {
  analysis_dir <- paste(paste0(output_dir_location,'/analysis_output/'))
  figure_dir <- paste(analysis_dir, Sys.Date(), sep="")
if (file.exists(output_dir_location)){print(paste0('generating folders in ',output_dir_location))} else {dir.create(file.path(output_dir_location))
  print('generating output directory')
}
  if (file.exists(analysis_dir)){print(paste(analysis_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(analysis_dir))
  print('generating analysis directory')
}
if (file.exists(figure_dir)){print(paste(figure_dir, 'exists, overwriting files if present'))} else {dir.create(file.path(figure_dir))
  print('generating figure directory')
}
return(figure_dir)
}

#' load seq2science countable
#'
#' This function loads the seq2science counts.tsv file
#' @param seq2science_dir Path to the directory where the seq2science results are located
#' @return counts in a dataframe
#' @export
load_s2s_counts_data <- function(seq2science_dir) {
  counts_folder <- paste0(seq2science_dir, 'results/counts/')
  if (file.exists(counts_folder)){
    counts_file <- list.files(counts_folder,pattern = 'counts.tsv')
    if (rlang::is_empty(counts_file)){print(paste0('cant find a file ending with count.tsv within ',counts_folder))}else{
      counts_file_path <- paste0(counts_folder,counts_file)
      print(paste0('loading the countmatrix ',counts_file_path))
      count_matrix <- read.table(file=counts_file_path, 
                                 row.names = 1, 
                                 header = T,
                                 sep = '\t')
      return(count_matrix)
    }
  } else {print(paste0('cant find ', counts_folder ))}
}
#' load seq2science sample file
#'
#' This function loads the seq2science samples.tsv file
#' @param seq2science_dir Path to the directory where the seq2science results are located
#' @return sample file as a dataframe
#' @export
load_s2s_sample_file <- function(seq2science_dir) {
  log_folder <- paste0(seq2science_dir, 'results/log/')
  if (file.exists(log_folder)){
    sample_file <- list.files(log_folder,pattern = 'samples.tsv')
    if (rlang::is_empty(sample_file)){print(paste0('cant find a file ending with samples.tsv within ',log_folder))}else{
      sample_file_path <- paste0(log_folder,sample_file)
      print(paste0('loading the samplefile ',sample_file_path))
      sample_file <- read.csv(file=sample_file_path, sep = '\t')
      return(sample_file)
    }
  } else {print(paste0('cant find ', log_folder ))}
}