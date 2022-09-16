     # Startup R in Ubuntu. Note, to escape out of a running program on Linux use < Ctrl-c >.
     R
     
     # < Ctrl-z > jumps out into the bash shell; < ps > in the shell shows processes; < fg > in the shell goes back into R 
     # (One could start another R process and jump between them, Google for more info.)
     
     # Now in R (see below for GitHub install)
     options(width = 140) # Default command line width in Linux-alikes is too short - adjust 140 to fit your current window size
     install.packages(c('sys', 'askpass', 'jsonlite', 'mime', 'openssl', 'R6', 'curl', 'httr', 'remotes', 'TMB'))
     
     # Installing TMB from GitHub also works (or for any other package that is on GitHUb)
     # There is less scrolling of C++ code on the screen with the GitHub install.
     # install.packages(c('sys', 'askpass', 'jsonlite', 'mime', 'openssl', 'R6', 'curl', 'httr', 'remotes'))
     # Sys.setenv(GITHUB_PAT='< Put your own GITHUB_PAT here >')  #  May work without this, at least for awhile.
     # remotes::install_github("kaskr/adcomp/TMB", INSTALL_opts = "--no-staged-install")
     
     # Test TMB
     
     library(TMB)
     runExample('simple')
     opt
     
     # Create an out-of-bounds error and verify the gdbsource() gives the correct line where the error occurs
     # Compiling 'simpleError.cpp' works:
     if(file.exists('simpleError.o')) file.remove(c('simpleError.o')); if(file.exists('simpleError.so')) file.remove(c('simpleError.so'))
     compile('simpleError.cpp', "-O0 -g")
     
     # But sourcing 'simpleError.R' will crash R
     source('simpleError.R')
     
     # Therefore, get back into R and use gdbsource() in Linux to find the line of code (30) that has the error
     library(TMB)
     gdbsource('simpleError.R')
   