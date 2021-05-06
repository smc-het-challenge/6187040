FROM ubuntu

# install R
FROM r-base:latest

# install R packages
RUN echo 'install.packages(c("data.table","gtools"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \
    && Rscript /tmp/packages.R

# install python
RUN apt-get update
RUN apt-get install -y python python-dev python-distribute python-pip

# install python packages contain in requirements
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app
COPY requirements.txt /usr/src/app/
RUN pip install -r requirements.txt
COPY . /usr/src/app

# add our own scripts 
WORKDIR /opt/QS
ADD montecnZaI5.R /opt/QS/montecnZaI5.R
ADD merge_vcfBat.py /opt/QS/merge_vcfBat.py
