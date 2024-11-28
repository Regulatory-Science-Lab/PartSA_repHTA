# Script to install required packages - run if these packages are not currently
# isntalled

install.packages(c('pacman', 'devtools'))
library('devtools')
install_github('DARTH-git/dampack') 