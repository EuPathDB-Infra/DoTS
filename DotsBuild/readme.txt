The DotsBuild component is used to build the DoTS transcript index.

INSTALLATION

First Time
1. Create a directory which will house the components from cvs, eg:
     % mkdir $HOME/components

2. Set your COMPONENTS_HOME environment variable to that location (you should
   add this to your .login or .post.cshrc file):
     % setenv COMPONENTS_HOME $HOME/components

3. Create the bin/ directory where the executables will reside.  This 
   directory should be considered volatile, ie, only housing things that
   are installed and can be deleted by an uninstall.  It should not hold
   hand-edited files. 
     % mkdir $HOME/bin

4. Add this directory to the front of your PATH (you should
   add this to your .login or .post.cshrc file):
     % setenv PATH $HOME/bin:$PATH

5. Create the lib/perl/ directory where the perl libraries will reside. This 
   directory should be considered volatile, ie, only housing things that
   are installed and can be deleted by an uninstall.  It should not hold
   hand-edited files. 
     % mkdir -p $HOME/lib/perl

6. Add that directory to the front of your PERL5LIB (you should
   add this to your .login or .post.cshrc file):
     % setenv PERL5LIB $HOME/lib/perl:$PERL5LIB

7. Set up SSH 2 (if its not already set up on local machine and alpha)
  - on the local machine:
     % cd ~/.ssh
     % ssh-keygen -t dsa
     % scp id_dsa.pub $USER@alpha.genomics.liniac.upenn.edu:~/.ssh
  - on alpha:
     % cd ~/.ssh
     % cat id_dsa.pub > authorized_keys2

8. Setup GUS config file, eg $HOME/.gus.cfg (see sample file below)

9. Set GUS_CFG to point to the config file (add this to your .login):
     % setenv GUS_CFG $HOME/.gus.cfg
    

Each Time
1. checkout the DotsBuild install script
     % cd $COMPONENTS_HOME
     % cvs co DotsBuild/install

2. Run the install script
     % cd $COMPONENTS_HOME
     % DotsBuild/install $HOME/bin $HOME/lib/perl $HOME/test -checkoutfirst

3. Fire off ssh agent (if not already running)
    % eval `ssh-agent`
    % ssh-add

Installing on the Liniac
1. set up PATH and PERL5LIB on alpha just as you did locally

2. copy over the installation:
     % cd $HOME
     % scp bin lib/perltest $USER@alpha.genomics.liniac.upenn.edu:$HOME



## 
Example config file for GUSdev
##

gus/database            GUSdev
gus/login               <db_login>
gus/password            <db_pasword>
gus/user                brunkb
gus/group               CBIL
gus/project             GUS