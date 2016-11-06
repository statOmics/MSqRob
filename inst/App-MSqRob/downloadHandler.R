  output$downloadData <- downloadHandler(filename = function() { paste0(input$project_name,".zip") }, #default name
                                         content = function(file){

                                           #Works, but all folders on top are also included in the zip:

                                           #!!!Temporary folder will be removed, careful when changing this!!!
                                           temppath <- file.path(getwd(), "temp")
                                           dir.create(temppath)

                                           proteins <- outputlist()$RData$proteins
                                           models <- outputlist()$RData$models
                                           save(proteins,models, file=file.path(temppath, "models.RData"))

                                           results <- outputlist()$results
                                           
                                           if(Sys.info()['sysname']=="Windows"){
                                             names_res_xlsx <- character()
                                             for(i in 1:length(results)){
                                               names_res_xlsx[i] <- paste0("results",i,".xlsx")
                                           xlsx::write.xlsx(results[[i]], file = file.path(temppath, names_res_xlsx[i]), col.names = TRUE, row.names = TRUE)
                                             }
                                             files <- file.path(temppath, "models.RData")
                                             for(i in 1:length(results)){
                                               files[(i+1)] <- file.path(temppath, names_res_xlsx[i])
                                             }
                                             zip(zipfile=file, files=files)
                                           }else{
                                           openxlsx::write.xlsx(results, file = file.path(temppath, "results.xlsx"), colNames = TRUE, rowNames = TRUE)
                                            zip(zipfile=file, files=c(file.path(temppath, "models.RData"),file.path(temppath, "results.xlsx")))
                                             }
                                           
 
                                           #!!!Removes temporary folder, careful when changing this!!!!
                                           unlink(temppath, recursive = TRUE)

                                         },
                                         contentType = "application/zip"
  )

