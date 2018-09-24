centre <- function(x, type) {
   switch(type,
            mean = mean(x),
            median = median(x),
            trimmed = mean(x, trim = .1))
  }
 x <- rcauchy(10)
 centre(x, "mean")

 centre(x, "median")
 centre(x, "trimmed")

covchoose <- function(covmodel) {
  switch(covmodel,
       exponential )
}
?eval


get.ss <- function(x) {
	return(sum(x ^ 2))
}
get.ss.sqrt <- function(x) {
	return(sqrt(sum(x ^ 2)))
}

get.something <- function(x, type) {
	switch(type,
		ss = eval(get.ss(x)),
		sqrtss = eval(get.ss.sqrt(x)))
}
get.something(x, type = "sqrtss")

x <- 1:10
get.ss(x)
eval(get.ss(x))
