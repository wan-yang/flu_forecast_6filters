## Script for initializing some parameters for the functions
## function to asign data (obs), find the first and last dates of a season
## The dates are matched with the CDC Flu View website
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

Fn_dates=function(season){
  if(season=="03-04"){
    weeks=1:53; # 53 weeks
    start_date=as.Date("2003/10/4");
    end_date=as.Date("2004/10/2"); 
  } else if (season=="04-05"){
    weeks=54:105;# October 3, 2004 onwards
    start_date=as.Date("2004/10/9"); # end of the week
    end_date=as.Date("2005/10/1");
  } else if (season=="05-06"){
    weeks=106:157  # October 3, 2005 onwards 
    start_date=as.Date("2005/10/8"); # end of the week
    end_date=as.Date("2006/09/30");
  } else if (season=="06-07"){
    weeks=158:209   # October 2, 2006
    start_date=as.Date("2006/10/7"); # end of the week
    end_date=as.Date("2007/9/29");
  } else if (season=="07-08"){
    weeks=210:261;
    start_date=as.Date("2007/10/6"); # end of the week
    end_date=as.Date("2008/9/27");
  } else if (season=="08-09") {
    weeks=262:291; # 30 Weeks
    start_date=as.Date("2008/10/4"); # end of the week
    end_date=as.Date("2009/4/25"); # week 16, end of the week, before the pandemic
  } else if (season=="09-10"){
    weeks=315:366;  
    start_date=as.Date("2009/10/10"); # end of the week
    end_date=as.Date("2010/10/2");
  } else if (season=="10-11"){
    weeks=367:418;
    start_date=as.Date("2010/10/9"); # end of the week
    end_date=as.Date("2011/10/1");
  } else if (season=="11-12"){
    weeks=419:470;
    start_date=as.Date("2011/10/8"); # end of the week
    end_date=as.Date("2012/9/29");  
  } else {
    weeks=NA;
    start_date=NA; # end of the week
    end_date=NA;  
  }
  rec=list(weeks=weeks,start_date=start_date,end_date=end_date)
  rec;
}

# function to calculation the mode
MODE <- function(x) {
  md=count(x)
  mode=md$x[which.max(md$freq)];
  count=max(md$freq);
  c(mode,count)
}