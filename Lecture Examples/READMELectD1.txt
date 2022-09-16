Instructions to install the TMB R package (https://github.com/kaskr/adcomp) running under command-line Ubuntu using WSL (Windows Subsystem for Linux).
This installation allows TMB::gdbsource() to properly debug Cpp files while using a PC running Windows 10 as the main OS.
It also provides R on an Ubuntu installation with GDB debugging software for any other purpose. For those with a different Ubuntu installation such as a dual-boot PC or VirtualBox, jumping past the WSL instuctions will give a guide to installing the folowing apps on Ubuntu release 18.04: curl, libssl-dev (for the httr R package), GDB, and the most recent version of R. 

WSL saves all the software in a single folder which can be used for computational reproducibility, exporting to another machine, or sharing with a colleague. The path to the WSL folder is given below. See also:

https://www.hanselman.com/blog/easily-move-wsl-distributions-between-windows-10-machines-with-import-and-export

     # You will need to be an administrator on your Windows machine to install WSL.
     
     # Run PowerShell as Administrator (use search (the magnifier glass icon) and type < power > and click on 'Run as Administrator'
     # Right clicking in PowerShell will insert what previously has been copied into the clipboard inside or outside or the PowerShell.
     
     # If the PowerShell starts in C:\WINDOWS\system32 the change to c:\, if desired
     cd ../.. 
     
     # In the PowerShell first enable wsl to be installed in Windows:
     PS> Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
     
     # For help with the wsl command do:
     PS> wsl --help
     
     # Get a list of valid distributions. Ubuntu 20.04 is available but may have issues.
     PS> wsl  --list --online
     
     # Install Ubuntu-18.04 into wsl by entering  :
     PS> wsl --install -d Ubuntu-18.04 
     
     # A new window will open. Enter a user name and password. A (very) short password be accepted (a single space will do) 
     #   and will be faster to type when needing elevated rights.
     # Exit out of the new window by typing < exit > and then start WSL in the PowerShell. < exit > also exits out of WSL.
     PS> wsl
     
     # If nothing happens when typing 'wsl' or the wrong distribution is opened, enter:
     PS> wsl -l
     PS> wsl -s Ubuntu-18.04
     
     # < wsl -s > sets the given distribution as the default. Now try WSL again:
     PS> wsl 
     
     # Check Ubuntu version (note that this flavor of Ubuntu is 'Bionic')
     lsb_release -a
     
     # In Windows File Explorer create the folder "TMB_Debug" on the C: drive. Enter the folder and leave File Explorer open.
     # Copy the 'simpleError.cpp' and simpleError.R files from the R_and_Cpp folder in this repo to C:\TMB_Debug
     
     # All the Windows' drives will be mounted in /mnt on Ubuntu, e.g. 'c' for 'C:' drive:
     ls ../../mnt
     
     # Create a change directory command in .bashrc to run at startup. '.bashrc' needs to be in the home folder (~).
     cd ~
     echo 'cd /mnt/c/TMB_Debug' > .bashrc
     ls -al
     
     # For now, change to this directory manually
     cd /mnt/c/TMB_Debug
     
     # Before installing the apps, including R, here is some additional information.
       
     # WSL files are in an \AppData\Local\Packages folder by user, e.g. mine are here:
     C:\Users\John\AppData\Local\Packages
     
     # The name (sort by most recently 'Date modified') will appear similar to:
     CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc
     
     # This folder can be renamed, backed up, saved for computational reproducibility, and exported to a colleague.  
     # (If you also use Docker, you may have exit out of Docker before working with the folder.)
     https://www.hanselman.com/blog/easily-move-wsl-distributions-between-windows-10-machines-with-import-and-export
     
     # WSL can also be opened by using the Windows' search feature: type < wsl > and click on 'Run as Administrator'.  
     # (There will be a different default color scheme.)
     
     
     # Make yourself root to do the following installs using the password set above.
     # The '$' prompt will change to '#', type < exit > to close. (Mnemonic: 'The pound is stronger than the dollar.')
     sudo su
     
     # First do an update
     apt-get update
     
     # Install 'curl' (The handy '-y' option auto answers 'Yes' to prompts to continue or stop.)
     apt-get -y install curl
     apt-get -y install libcurl4-openssl-dev
          
     # Install libssl-dev for the R package 'httr'
     apt-get -y install libssl-dev
          
     # Install gdb following: https://stackoverflow.com/questions/10255082/installing-r-from-cran-ubuntu-repository-no-public-key-error
     # The first command may take some time or perhaps will fail for a time - continue to retry - install R and return here if needed
     gpg --keyserver pgp.mit.edu --recv-key 51716619E084DAB9
     
     gpg -a --export 51716619E084DAB9 > jranke_cran.asc 
     apt-key add jranke_cran.asc 
     apt -y install gdb
     
     # Install BLAS, LAPACK, and FORTRAN packages
     apt-get -y install libblas-dev liblapack-dev
     apt-get -y install gfortran
     
     
     # CRAN instructions for installing R in Ubuntu: https://cran.r-project.org/  
     # On the CRAN website select 'Ubuntu' on the line: < Download R for Linux (Debian, Fedora/Redhat, Ubuntu) >
     # < install build-essential> is not in the CRAN instuctions but was needed by me and others on the internet.
     apt update -qq
     apt -y install build-essential
     
     # install two helper packages we need
     apt -y install --no-install-recommends software-properties-common dirmngr
     
     # add the signing key (by Michael Rutter) for these repos
     # To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
     # (The verify key command doesn't work for me, but isn't needed to continue.)
     # Fingerprint: 298A3A825C0D65DFD57CBB651716619E084DAB9
     wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
     
     # add the R 4.0 repo from CRAN
     # Here we use lsb_release -cs to access which Ubuntu flavor you run: one of “impish”, “hirsute”, “focal”, “bionic”, …
     add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
          
     # Then run
     apt -y install --no-install-recommends r-base
     
   