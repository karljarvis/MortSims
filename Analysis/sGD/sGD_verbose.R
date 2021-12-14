sGD_verbose = function (genind_obj, xy, dist.mat, radius, min_N, NS_ans = F, 
                        GD_ans = T, NHmat_ans = F, genout_ans = F, file_name = NULL, 
                        NeEstimator_dir = NULL) 
{
  if (packageVersion("adegenet") < "2.0.0") {
    stop("Please install the latest version of the adegenet package (>= 2.0.0)")
  }
  if (packageVersion("hierfstat") < "0.04.15") {
    stop("Please install the development version of the hierfstat package (>= 0.04.15) \n       First, install the devtools package, and then run:\n               library(devtools)\n               install_github(\"jgx65/hierfstat\")\n       After installing hierfstat, please restart R before running sGD.")
  }
  if (NS_ans == F & GD_ans == F) {
    stop("At least one of the following must be TRUE: NS_ans or GD_ans")
  }
  cat("Reading input files...\n")
  df_dat = genind2df(genind_obj)
  allele.digits = nchar(genind_obj@all.names[[1]][1])
  numloci = length(names(genind_obj@all.names))
  numindivs = dim(genind_obj@tab)[1]
  OS = as.character(Sys.info()["sysname"])
  NH_genepop_filename = paste(file_name, "_genepop.gen", sep = "")
  if (is.null(file_name) == FALSE) {
    NH_summary_filename = paste(file_name, "_sGD.csv", sep = "")
  }
  if (is.numeric(min_N) == F) {
    stop("min_N must be an integer")
  }
  if (is.numeric(xy[, 2]) == F) {
    stop("The second column of the xy file must be a numeric x coordinate (e.g. longitude)")
  }
  if (is.numeric(xy[, 3]) == F) {
    stop("The third column of the xy file must be a numeric y coordinate (e.g. latitude)")
  }
  if (nrow(xy) != numindivs) {
    stop("The number of rows in the xy file does not match the number of individuals in the genepop file")
  }
  if (nrow(dist.mat) != numindivs) {
    stop("The number of rows in the dist.mat does not match the number of individuals in the genepop file")
  }
  if (ncol(dist.mat) != numindivs) {
    stop("The number of columns in the dist.mat does not match the number of individuals in the genepop file")
  }
  if (is.null(NeEstimator_dir) == F) {
    if (OS == "Windows") {
      if (file.exists(file.path(NeEstimator_dir, "Ne2.exe")) == 
            F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
    }
    if (OS == "Linux") {
      if (file.exists(file.path(NeEstimator_dir, "Ne2L")) == 
            F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
      if (file.exists(file.path(NeEstimator_dir, "NeEstimator.jar")) == 
            F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
    }
    if (OS == "Darwin") {
      if (file.exists(file.path(NeEstimator_dir, "Ne2M")) == 
            F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
      if (file.exists(file.path(NeEstimator_dir, "NeEstimator.jar")) == 
            F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
    }
  }
  if (is.null(NeEstimator_dir) == T & NS_ans == T) {
    stop("NS_ans is TRUE, however, you have not specified the location of NeEstimator_dir")
  }
  cat("Input summary:\n")
  cat(paste("\t individuals:", numindivs, "\n"))
  cat(paste("\t loci:", numloci, "\n"))
  cat(paste("\t neighborhood radius:", radius, "\n"))
  cat(paste("\t minimum sample size:", min_N, "\n"))
  if (is.null(file_name) == FALSE) {
    cat(paste("\t output file:", NH_summary_filename, "in", 
              getwd(), "\n"))
  }
  header_length = 1 + numloci
  genepop_header = "Genepop_file"
  write.table(genepop_header, NH_genepop_filename, sep = "\t", 
              quote = F, col.names = F, row.names = F)
  for (locus in 1:numloci) {
    write.table(names(genind_obj@all.names)[locus], NH_genepop_filename, 
                sep = "\t", quote = F, col.names = F, row.names = F, 
                append = T)
  }
  cat("Determining neighborhood membership from dist.mat and radius...\n")
  neighborhood_mat = ifelse(dist.mat < radius, 1, 0)
  neighborhood_N = colSums(neighborhood_mat)
  valid_pops = as.numeric(which(neighborhood_N >= min_N))
  npops = length(valid_pops)
  for (pop in valid_pops) {
    write.table("POP", NH_genepop_filename, sep = "\t", quote = F, 
                col.names = F, row.names = F, append = T)
    output = df_dat[which(neighborhood_mat[pop, ] == 1), 
                    ]
    indivIDs = paste("P", pop, "I", c(1:nrow(output)), sep = "")
    output = data.frame(indivIDs, ",", output[2:ncol(output)])
    write.table(output, NH_genepop_filename, sep = "\t", 
                quote = F, col.names = F, row.names = F, append = T)
  }
  NH_summary = data.frame(dimnames(genind_obj@tab)[[1]], xy[, 
                                                            c(2, 3)], neighborhood_N, stringsAsFactors = F)
  names(NH_summary) = c("Indiv_ID", "X", "Y", "N")
  NH_summary$index = c(1:numindivs)
  if (GD_ans == T) {
    cat("Calculating genetic diversity indices for neighborhoods...\n")
    NH_genind = read.genepop(NH_genepop_filename, ncode = allele.digits, 
                             quiet = T)
    NH_Ar = round(colMeans(allelic.richness(NH_hierfstat)$Ar), 
                  4)
    NH_stats = basic.stats(NH_hierfstat, digits = 4)
    NH_Ho = round(colMeans(NH_stats$Ho), 4)
    NH_Hs = round(colMeans(NH_stats$Hs), 4)
    NH_FIS = round(1 - NH_Ho/NH_Hs, 4)
    GD_output = data.frame(NH_summary$index[valid_pops], 
                           NH_Ar, NH_Hs, NH_Ho, NH_FIS, stringsAsFactors = F)
    names(GD_output) = c("index", "Ar", "Hs", "Ho", "FIS")
    NH_summary = merge(NH_summary, GD_output, by = "index", 
                       all = T)
    print(numindivs)
  }
  if (NS_ans == T) {
    write.table(1, "Ne2_input.txt", sep = "\t", col.names = F, 
                row.names = F, quote = F)
    write.table(3, "Ne2_input.txt", sep = "\t", col.names = F, 
                row.names = F, quote = F, append = T)
    write.table(cbind(0.1, 0.05, 0.02), "Ne2_input.txt", 
                sep = "\t", col.names = F, row.names = F, quote = F, 
                append = T)
    write.table(NH_genepop_filename, "Ne2_input.txt", sep = "\t", 
                col.names = F, row.names = F, quote = F, append = T)
    write.table("Ne2_output.txt", "Ne2_input.txt", sep = "\t", 
                col.names = F, row.names = F, quote = F, append = T)
    cat("Calculating NS for neighborhoods...\n")
    if (OS == "Windows") {
      file.copy(file.path(NeEstimator_dir, "Ne2.exe"), 
                getwd())
      system("Ne2.exe  m:Ne2_input.txt", show.output.on.console = F, 
             ignore.stdout = T, ignore.stderr = T)
      file.remove("Ne2.exe")
    }
    if (OS == "Darwin") {
      file.copy(file.path(NeEstimator_dir, "Ne2M"), getwd())
      file.copy(file.path(NeEstimator_dir, "NeEstimator.jar"), 
                getwd())
      system("./Ne2M  m:Ne2_input.txt", ignore.stdout = T, 
             ignore.stderr = T)
      file.remove("Ne2M")
      file.remove("NeEstimator.jar")
    }
    if (OS == "Linux") {
      file.copy(file.path(NeEstimator_dir, "Ne2L"), getwd())
      file.copy(file.path(NeEstimator_dir, "NeEstimator.jar"), 
                getwd())
      system("./Ne2L  m:Ne2_input.txt", ignore.stdout = T, 
             ignore.stderr = T)
      file.remove("Ne2L")
      file.remove("NeEstimator.jar")
    }
    LDNe_output = readLines("Ne2_output.txt")
    LDNe_datalines = grep("Estimated Ne", LDNe_output)
    LDNe_data = unlist(strsplit(LDNe_output[LDNe_datalines], 
                                "\\s+"))
    LDNe_estimates = data.frame(matrix(LDNe_data, nrow = npops, 
                                       ncol = 7, byrow = T)[, 4:7], stringsAsFactors = F)
    LDNe_estimates = data.frame(NH_summary$index[valid_pops], 
                                LDNe_estimates, stringsAsFactors = F)
    names(LDNe_estimates) = c("index", "NS_ex0.10", "NS_ex0.05", 
                              "NS_ex0.02", "NS_ex0.00")
    NH_summary = merge(NH_summary, LDNe_estimates, by = "index", 
                       all = T)
    file.remove("Ne2_input.txt")
    file.remove("Ne2_output.txt")
    if (genout_ans == FALSE) {
      file.remove(NH_genepop_filename)
    }
  }
  for (col in c(3:ncol(NH_summary))) {
    NH_summary[, col] = suppressWarnings(as.numeric(NH_summary[, 
                                                               col]))
  }
  if (is.null(file_name) == FALSE) {
    cat("Appending results to neighborhood summary file...\n")
    write.table(NH_summary, NH_summary_filename, row.names = F, 
                sep = ",", na = "")
  }
  if (NHmat_ans == TRUE) {
    cat("Writing neighborhood membership matrix to file...\n")
    write.table(neighborhood_mat, paste(file_name, "_neighborhood_mat.csv", 
                                        sep = ""), row.names = F, col.names = F, sep = ",", 
                na = "")
  }
  return(NH_summary)
  cat("Processing complete.\n")
}