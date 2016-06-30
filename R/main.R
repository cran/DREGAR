dregar <- function  ( data,  da = 0,  ar = 0 ,  
                      mselection  =  4       ,
                      normalize = FALSE      ,  
                      penalized = TRUE      ,  
                      iteration = 15) 
{ 
  result = regarma(x = data[,-c(1,2)],y = data[,2],
                   ar = da , ma = ar,
                   normalize = normalize , debug = 0,
                   rep = iteration , pen = penalized,
                   criteria = mselection 
  )
  
  return(result)
} 


