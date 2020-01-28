

# Function to import the data ---------------------------------------------

generate_dataset <- function(data_folder = 'NKI_Rockland/',
                             y_filename = 'clinical_information.txt',
                             y_column = 'WASI_FULL_4',
                             output_filename = 'NKIdata.RData',
                             output_folder = ".",
                             ext_save = FALSE){

  # Response variable
  y_df <- read.csv(y_filename, colClasses = c('Subject' = 'character'))[c('Subject', y_column)]

  # Exclude obs with no response
  y_df <- na.exclude(y_df)

  # Ids for the patients who have the response
  id_y <- y_df$Subject

  # Function to retrieve the ids for the patients who have a certain 'pattern'
  id_pattern <- function(pattern){
    files <- list.files(data_folder, pattern = pattern)
    sapply(files, FUN = function(file){
      id <- strsplit(file, split = '_', fixed = TRUE)[[1]][1]
    })
  }

  # Ids for the patients who have DTI
  id_dti <- id_pattern('_DTI_connectivity_matrix_file.txt')

  # Ids for the patients who have noGSR fcMRI
  id_fcmri_nogsr <- id_pattern('_fcMRI_noGSR_connectivity_matrix_file.txt')

  # Ids for the patients who have response, DTI and noGSR fcMRI
  id_final <- Reduce(intersect, list(id_y, id_dti, id_fcmri_nogsr))

  # Dataset initialization
  data_out <- list(y = y_df$WASI_FULL_4)

  # Retrieve all the files' names for the DTI matrices
  dti_filenames <- sapply(id_final, function(id){
    paste0(data_folder, id, '_DTI_connectivity_matrix_file.txt')
  })

  # Retrieve all the files' names for the noGSR fcMRI matrices
  fcmri_filenames <- sapply(id_final, function(id){
    paste0(data_folder, id, '_fcMRI_noGSR_connectivity_matrix_file.txt')
  })

  # Import the DTI matrices
  data_out$structural <- lapply(dti_filenames, function(d) as.matrix(read.table(d)))

  # Import the noGSR fcMRI matrices
  data_out$functional <- lapply(fcmri_filenames, function(f) as.matrix(read.table(f)))

  # Optional saving
  if(ext_save){
    save(data_out, file = paste0(output_folder, output_filename))
  }

  return(data_out)

}



# Data import -------------------------------------------------------------

nki <- generate_dataset(data_folder = 'NKI_Rockland/',
                        y_filename = 'clinical_information.txt',
                        y_column = 'WASI_FULL_4',
                        output_filename = 'NKIdata.RData',
                        output_folder = ".",
                        ext_save = FALSE)



# Label comparison --------------------------------------------------------

### DTI ###

# DTI_region_names_full_file comparison
dti_full_files <- list.files(pattern='_DTI_region_names_full_file.txt', full.names = FALSE)
comp_full_files <- sapply(dti_full_files, FUN = function(dti) {
  all.equal(readLines('1013090_DTI_region_names_full_file.txt'), readLines(dti))
}
)
comp_full_files
#many differences

# DTI_region_names_abbrev_file comparison
dti_abbrev_files <- list.files(pattern='_DTI_region_names_abbrev_file.txt', full.names = FALSE)
comp_abbrev_files <- sapply(dti_abbrev_files, FUN = function(dti) {
  all.equal(readLines('1013090_DTI_region_names_abbrev_file.txt'), readLines(dti))
}
)
all(comp_abbrev_files)
#all equal

# DTI_region_xyz_centers_file comparison
dti_xyz_files <- list.files(pattern='_DTI_region_xyz_centers_file.txt', full.names = FALSE)
comp_xyz_files <- sapply(dti_xyz_files, FUN = function(dti) {
  all.equal(readLines('1013090_DTI_region_xyz_centers_file.txt'), readLines(dti))
}
)
all(comp_xyz_files)
#all equal


### fRMI (noGSR) ###

# fcRMI_region_names_full_file comparison
fcmri_full_files <- list.files(pattern='_fcMRI_noGSR_region_names_full_file.txt', full.names = FALSE)
comp_full_files <- sapply(fcmri_full_files, FUN = function(fcmri) {
  all.equal(readLines('1013090_fcMRI_noGSR_region_names_full_file.txt'), readLines(fcmri))
}
)
comp_full_files
#many differences

# fcMRI_noGSR_region_names_abbrev_file comparison
fcmri_abbrev_files <- list.files(pattern='_fcMRI_noGSR_region_names_abbrev_file.txt', full.names = FALSE)
comp_abbrev_files <- sapply(fcmri_abbrev_files, FUN = function(fcmri) {
  all.equal(readLines('1013090_fcMRI_noGSR_region_names_abbrev_file.txt'), readLines(fcmri))
}
)
all(comp_abbrev_files)
#all equal

# fcMRI_noGSR_region_xyz_centers_file comparison
fcmri_xyz_files <- list.files(pattern='_fcMRI_noGSR_region_xyz_centers_file.txt', full.names = FALSE)
comp_xyz_files <- sapply(fcmri_xyz_files, FUN = function(fcmri) {
  all.equal(readLines('1013090_fcMRI_noGSR_region_xyz_centers_file.txt'), readLines(fcmri))
}
)
all(comp_xyz_files)
#all equal
