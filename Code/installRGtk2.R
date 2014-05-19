windows_gtk_path <- function(){
    file.path(system.file(package = "RGtk2"), "gtk", .Platform$r_arch)
}
configure_gtk_theme <- function(theme) {
  ## Only applies to Windows so far
  config_path <- file.path(system.file(package = "RGtk2"), "gtk",
                           .Platform$r_arch, "etc", "gtk-2.0")
  dir.create(config_path, recursive = TRUE)
  writeLines(sprintf("gtk-theme-name = \"%s\"", theme),
             file.path(config_path, "gtkrc"))
}

windows32_config <- list(
    source = FALSE,
    gtk_url = "http://ftp.gnome.org/pub/gnome/binaries/win32/gtk+/2.22/gtk+-bundle_2.22.1-20101227_win32.zip",
    installer = function(path) {
        gtk_path <- windows_gtk_path()
        ## unzip does this, but we want to see any warnings
        dir.create(gtk_path, recursive = TRUE) 
        unzip(path, exdir = gtk_path)
        configure_gtk_theme("MS-Windows")})
        
windows64_config <- windows32_config
windows64_config$gtk_url <- "http://ftp.gnome.org/pub/gnome/binaries/win64/gtk+/2.22/gtk+-bundle_2.22.1-20101229_win64.zip"
           
         
if (.Platform$OS.type == "windows") {
    if (.Platform$r_arch == "i386"){
        config <- windows32_config
    } else {
        config <- windows64_config
    }
}

installer <- config$installer
dep_url <- config$gtk_url
dep_name <- "GTK+"
gtk_web <- "http://www.gtk.org"
# path <- file.path(tempdir(), basename(sub("\\?.*", "", dep_url))) 
ojtpath <- "C:/Users/jang/Documents/OJT"
path <- file.path(ojtpath, "R", "Programs", basename(sub("\\?.*", "", dep_url)))
installer(path)