FROM rocker/rstudio:4.3.2

ENV USER=rstudio
ENV PASSWORD=endo

COPY . /C:/Users/Erik/Desktop/Programming/R/Bio/Endo
WORKDIR /C:/Users/Erik/Desktop/Programming/R/Bio/Endo

EXPOSE 8787
