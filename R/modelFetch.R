
#' Load Pre-trained Models from Github
#' 
#' Use REST API of the github and download pre-trained models onto the 
#' local machine
#' 
#' @param url Main API url
#' @param path Repository path of the model 
#' @param version Version number
#' @param cancer.type Cancer type
#' @param token Authentication token; if `NA`, access attempted without one
#' @param verbose Print signature names being downloaded
#' @param min.M Minimum no. of mutations
#' @param min.attr Mininum attribution
#' @param verbose Verbosity level
#' @param model.path Model path where downloaded data are to be placed 
#' dir <- modelFetch()
#' print(dir)
#' 
#' @export
modelFetch <- function(url = 'https://api.github.com', 
                       path = '/repos/mskcc/DeepSig/contents/inst/extdata/dlsig/',
                       version = 'v0.95', cancer.type = 'breast', token = NA, 
                       verbose = TRUE, model.path){
  if(is.na(model.path))
    model.path <- tempdir()
  req <- httr2::request(url) |>
    httr2::req_url_path(paste0(path, version, '/', cancer.type, '/'))
  if(!is.na(token)) req <- req |> httr2::req_auth_bearer_token(token)
  resp <- req |> httr2::req_perform()
  if((resp |> httr2::resp_status_desc())!='OK') 
    stop(paste0('API query for ',cancer.type,' failed'))
  
  z <- httr2::resp_body_json(resp)
  if(verbose) cat('Querying github API for ',cancer.type,' models...\n',sep='')
  pb <- txtProgressBar(style=3)
  for(i in rev(seq_along(z))){
    zi <- z[[i]]
    setTxtProgressBar(pb, 1-(i/length(z)))
    if(zi$type=='file'){
      download.file(zi$download_url, paste0(model.path,'/',zi$name), quiet=TRUE)
      next()
    }
    dir <- paste0(model.path,'/',zi$name)
    if(dir.exists(dir)) system(paste0('rm -rf ',dir))
    dir.create(dir,showWarnings=FALSE)
    req2 <- httr2::request(url) |>
      httr2::req_url_path(paste0(path, version, '/', cancer.type,'/',zi$name))
    if(!is.na(token)) req2 <- req2 |> httr2::req_auth_bearer_token(token)
    resp2 <- req2|> httr2::req_perform()
    if((resp2 |> httr2::resp_status_desc())!='OK') stop(paste0('API query for ',cancer.type,': ', zi$name,' failed'))
    z2 <- httr2::resp_body_json(resp2)
    for(zj in z2){
      if(zj$type=='file'){
        download.file(zj$download_url, paste0(dir,'/',zj$name), quiet=TRUE)
        next()
      }
      dir2 <- paste0(dir, '/', zj$name)
      if(dir.exists(dir2)) system(paste0('rm -rf ',dir2))
      dir.create(dir2, showWarnings=FALSE)
      req3 <- httr2::request(url) |>
        httr2::req_url_path(paste0(path, version, '/', cancer.type, '/', zi$name, '/',zj$name))
      if(!is.na(token)) req3 <- req3 |> httr2::req_auth_bearer_token(token)
      resp3 <- req3 |> httr2::req_perform()
      if((resp3 |> httr2::resp_status_desc())!='OK')
        stop(paste0('API query for ',cancer.type,': ', zi$name,'/',zj$name,' failed'))
      z3 <- httr2::resp_body_json(resp3)
      for(zk in z3){
        if(zk$type=='file'){
          download.file(zk$download_url, paste0(dir2, '/', zk$name), quiet=TRUE)
          next()
        } else{
          stop(paste0('Error reading model ',cancer.type))
        }
      }
    }
  }
  close(pb)
  return(model.path)
}