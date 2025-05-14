# Author:		E. Reichenberger
# Date:			Winter 24/25
# Purpose: Combine benchmarking files into one, make pretty images.

lib_path <- '.'
project <- ''

args <- commandArgs(trailingOnly = TRUE)
print(length(args))

nargs = 2
if (length(args) != nargs) {
  stop(paste(nargs, 'arguments must be supplied.'), call. = FALSE)
} else if (length(args) == nargs) {
	lib_path = args[1]
	project = args[2] 
}

library(tidyverse, lib.loc=lib_path)
library(ggrepel, lib.loc=lib_path)

benchmarks <- 'benchmarks/'
project_benchmarks <- paste0(benchmarks, project, '/')
print(project_benchmarks)

list_of_files <- list.files(path=project_benchmarks, pattern = ".benchmark", full.names=TRUE)

combined_data <- do.call(rbind, lapply(list_of_files, function(file) {
  df <- read.table(file, header = TRUE)  
  df$benchmark_application <- tools::file_path_sans_ext(basename(file))  
  
  return(df)  
}))

combined_data <- combined_data[, c(2, 3, 4, 7, 8, 9, 10, 11)]

# Define + Make directories
# ----------------------------------------------------------------
benchmark_dir <- paste0('data/endpoints/', project, '/analysis/report/benchmarks/')
benchmark_table_dir <- paste0(benchmark_dir, 'table/')
benchmark_fig_dir <- paste0(benchmark_dir, 'figures/')

dir.create(benchmark_dir, showWarnings = FALSE)
dir.create(benchmark_table_dir, showWarnings = FALSE)
dir.create(benchmark_fig_dir, showWarnings = FALSE)
# ----------------------------------------------------------------

# Write Benchmarking Data to File
# ----------------------------------------------------------------
combined_data_file <- paste0(benchmark_table_dir, project, '_benchmarks.txt')
write.table(combined_data, file=combined_data_file, quote=FALSE, row.names=FALSE, sep='\t', append=FALSE)
# ----------------------------------------------------------------

# Create Images
# ----------------------------------------------------------------
# reorder df
combined_data <- combined_data[order(combined_data$h.m.s),]

# time vs memory
time_memory <- paste0(benchmark_fig_dir, project, '_time_memory.png')
png(filename=time_memory, width=2700,height=2000,res=300)
print(ggplot(combined_data, aes(h.m.s, max_rss, colour=benchmark_application)) + 
  geom_point(show.legend = FALSE) +  # Remove legend for points
  geom_text_repel(aes(label = benchmark_application), show.legend = FALSE) +
	theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5, hjust=.5), 
        axis.text.y = element_text(size = 12, hjust = 1, angle = 45),  
        axis.title.x = element_text(size = 16, margin = margin(t = 20, b = 10)),  
        axis.title.y = element_text(size = 16, margin = margin(l = 5))) +
  scale_y_continuous(breaks = pretty(combined_data$max_rss, n = 15))
)
dev.off()

# memory read vs write
in_out_file <- paste0(benchmark_fig_dir, project, '_in_out.png')
png(filename=in_out_file, width=2700,height=2000,res=300)
print(ggplot(combined_data, aes(io_in, io_out, colour=benchmark_application)) + 
  geom_point() + 
  geom_text_repel(aes(label = benchmark_application), show.legend = FALSE, max.overlaps = 10) +
	theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5, hjust=.5), 
        axis.text.y = element_text(size = 12, hjust = 1, angle = 45),  
        axis.title.x = element_text(size = 16, margin = margin(t = 20, b = 10)),  
        axis.title.y = element_text(size = 16, margin = margin(l = 5))) + 
        scale_y_continuous(breaks = pretty(combined_data$io_out, n = 15)) +
        scale_x_continuous(breaks = pretty(combined_data$io_in, n = 10))
)
dev.off()
# ----------------------------------------------------------------
