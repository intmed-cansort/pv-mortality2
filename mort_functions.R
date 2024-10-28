# functions used for each cancer type
get_data<-function(cancer_type) {
  mort<-  haven::read_dta("data/mort.dta") |>    
    dplyr::filter(cancer_site1==cancer_type) |>  
    haven::as_factor() |>
    dplyr::mutate(stage = relevel(as.factor(stage), 3)) |>      # set base stage II
    dplyr::mutate(agedec = agedec-mean(agedec)) |>
    dplyr::mutate(pvgroup = relevel(as.factor(pvgroup), 5)) |>    # set base none
    dplyr::mutate(grade = relevel(as.factor(grade),2))  |>
    dplyr::mutate(race = relevel(race,4)) |>
    dplyr::mutate(rurality = relevel(rurality,1)) |>
    dplyr::mutate(poverty = relevel(poverty,3)) 
  return(mort)
}

get_fixed<-function(data,covar) {
  outcome<-"survival::Surv(survtime, event)"
  f<-as.formula(paste(outcome,paste(covar,collapse="+"),sep="~"))
  m_fixed<-survival::coxph(f, data, ties = "breslow")
  return(m_fixed)
}

get_random<-function(data,can_type,covar) {
  library(coxme)
  outcome<-"survival::Surv(survtime, event)"
  f<-as.formula(paste(outcome,paste(paste(covar[-length(covar)],collapse="+"),
                                    "cluster(pvgroup)", sep="+"),sep="~"))
  m_emfrail<-frailtyEM::emfrail(f, data)
  save_file<-paste0("data/","m_emfrail",can_type,".Rdata")
  save(m_emfrail,file=save_file)
  return(m_emfrail)
}

process_fixed<-function(data=m_fixed) {
  library(tidyverse)
  # Obtain center estimates and calculate standard error
  pvgroup <- grep("pvgroup", names(m_fixed$coefficients))
  nc <- length(pvgroup)
  var_pvgroup <- m_fixed$var[pvgroup, pvgroup]
  varmean <- sum(var_pvgroup)/ nc^2
  
  newse <- sqrt(c(varmean, diag(var_pvgroup) + varmean - 2/nc * apply(var_pvgroup, 1, sum)))
  # change contrast for hr's to average
  fe<-data.frame(beta = c(pvgroup5 = 0, m_fixed$coefficients[pvgroup]), sd = newse) |> 
    rownames_to_column("pvgroup") |> 
    mutate(pvgroup = substr(pvgroup,8,25)) |>
    mutate( pvgroup = base::replace(pvgroup, pvgroup=="5", "Negative"))  |>
    mutate(sbeta = sum(beta) / nc) |> 
    mutate(beta = beta - sbeta) |> 
    arrange(beta) |>
    mutate(ord = 1:n()) |> 
    mutate(pvgroup = as.character(pvgroup)) |> 
    mutate(ymin = beta - 1.96 * sd, ymax = beta + 1.96 * sd) |> 
    rename(id = pvgroup, estimate = beta) |> 
    select(-sd, -sbeta, -ord) |>
    mutate(type = "fixed effects") |> 
    mutate(ord = 1:n()) 
  return(fe)
}


process_random <- function(data=m_emfrail) {
  library(data.table)
  frailties<-as.data.table(frailtyEM:::summary.emfrail(data)[[12]])
  logz <- frailties[order(z),.(id,z,lower_q, upper_q)] |>
              dplyr::mutate_if(is.numeric, log) |> 
              dplyr::rename(estimate = z, ymin = lower_q, ymax = upper_q) |> 
              dplyr::mutate(type = "frailty") |>
              dplyr::mutate(ord = 1:dplyr::n())
return(logz)
}

combine_effects<-function(fe,logz) {
  ord_re <- data.frame(id = logz$id, ord = logz$ord)
  datt <- fe |>  
      right_join(ord_re,join_by(id==id)) |> 
      select(-ord.y) |> rename(ord=ord.x) |> 
      bind_rows(logz)
  return(datt)
}

table_fig <- function(tblname,size=1, tbl_title,type=can_type) {
  cancer_type<-c("breast cancer triple-negative","breast cancer her2+","breast cancer her2-,ER/PR+","colorectal","pancreas")
  gt::gt(dplyr::as_tibble(tblname,rownames="type")) |> 
    gt::tab_header(title = tbl_title,subtitle=cancer_type[type]) |>
    gt::gtsave(paste0("stats_",deparse(substitute(tblname)),
                      "_table",type,".png"), path = "img/")
  cowplot::ggdraw() + cowplot::draw_image(paste0("img/",
                                                 "stats_",deparse(substitute(tblname)),"_table",can_type,".png"), scale = size)
}

graph_effects_vert<-function(data=datt) {
  library(grid)
  library(data.table)
  datt_DT<-as.data.table(datt)
  datt_DT[,type :=(ifelse(type=="fixed effects","fe","re"))]
  
  dt<-dcast(datt_DT,
#            id + ord ~ type,
            id ~ type,
            value.var = c("estimate","ymin","ymax")) 
  dt$` ` <- paste(rep(" ", 35), collapse = " ")
  ests<-c("estimate_fe","estimate_re","ymin_fe","ymin_re","ymax_fe","ymax_re")
  dt[,(ests):=lapply(.SD, function(x) exp(x)),.SDcols = ests]
  
  #  form for estimates 1.2 (0.2, 2.9)\n1.6 (0.4, 3.5)
  dt$CI <- paste(sprintf("%.1f (%.1f, %.1f)", dt$estimate_fe, dt$ymin_fe, dt$ymax_fe),
                 sprintf("%.1f (%.1f, %.1f)", dt$estimate_re, dt$ymin_re, dt$ymax_re),
                 sep = "\n")
  setnames(dt, old = "id", new = "Gene")
  tm <- forestploter::forest_theme(base_size = 10,
                     refline_gp = gpar(lty="solid"),
                     ci_pch = c(15, 15),
                     ci_col = c("#377eb8", "#4daf4a"),
                     footnote_gp = gpar(col = "blue"),
                     legend_name = "Estimates",
                     legend_value = c("Fixed effects", "Random effects"),
                     legend_position = "none",
                     vertline_lty = c("dashed"),
                     vertline_col = c("#d6604d"),
                     # Table cell padding, width 4 and heights 3
                     core = list(padding = unit(c(4, 3), "mm")),
                     gp = gpar(fontface = "italic")
                     )
p<-forestploter::forest(dt[,c(1,8, 9)],
                          est = list(dt$estimate_fe,dt$estimate_re),
                          lower = list(dt$ymin_fe,dt$ymin_re), 
                          upper = list(dt$ymax_fe,dt$ymax_re),
                          ci_column = 2,
                          ref_line = 1, 
                          x_trans = "log",
                          ticks_at = c(0.25,0.5, 2,4),
                          arrow_lab = c("Lower risk", "Higher risk"),
                          nudge_y = 0.4,
                          xlim=c(.2,5),
                          theme=tm)
  return(p)
}
  
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
# https://stackoverflow.com/questions/23018256/printing-p-values-with-0-001
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  
  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('< %s', format(sig.limit)))
    if (x > .1)
      return(paste("=",roundr(x, digits = 2))) else
        return(paste("=",roundr(x, digits = digits)))
  }, sig.limit = sig.limit)
}


