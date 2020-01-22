basin_table <- function(basin, data, ...){
  station_rows <- map_df(unique(data$LOC),
                         ~filter(data, LOC == .x) %>% 
                           station_row()) %>% 
    arrange(order)
  bind_rows(basin_row(basin, data), station_rows)
}

basin_row <- function(basin, data, althyp = 1.1){
  n_LOC <- length(unique(data$LOC))
  if (n_LOC > 1) {
    basin_fit <- lm(log(value) ~ YEAR + LOC, data = data)
    basin_fit_sat <- lm(log(value) ~ YEAR * LOC, data = data)
    tibble(LOC = basin,
           intercept = common_intercept(basin_fit),
           slope = coef(basin_fit)["YEAR"],
           slope_se = vcov(basin_fit)["YEAR", "YEAR"] %>% sqrt(),
           relative_rate = (exp(slope) - 1) * 100,
           p_vs_sat = anova(basin_fit, basin_fit_sat)$"Pr(>F)"[2],
           power = (1 - pt(qt(.975,  basin_fit$df), ncp = log(althyp)/slope_se, df = basin_fit$df)) + pt(qt(.025,  basin_fit$df), ncp = log(althyp)/slope_se, df = basin_fit$df),
           limit = mean(data$limit),
           first_year = min(data$YEAR),
           n = nrow(data)) %>% 
      mutate(predicted_2018 = exp(intercept + slope * 2018),
             CI = paste0("(", round((exp(confint(basin_fit)[2, 1]) - 1) * 100, 2), ", ", round((exp(confint(basin_fit)[2, 2]) - 1) * 100, 2), ")"),
             ttt = ifelse(is.na(limit), NA, round((log(limit) - intercept) / slope) - 2018),
             ttt = ifelse(ttt < 0, Inf, ttt),
             n.all.lod = sum(data$all.lod))
  }
  else
  {
    basin_fit <- lm(log(value) ~ YEAR, data = data)
    tibble(LOC = basin,
           intercept = coef(basin_fit)["(Intercept)"],
           slope = coef(basin_fit)["YEAR"],
           slope_se = vcov(basin_fit)["YEAR", "YEAR"] %>% sqrt(),
           relative_rate = (exp(slope) - 1) * 100,
           p_vs_sat = NA,
           power = (1 - pt(qt(.975,  basin_fit$df), ncp = log(althyp)/slope_se, df = basin_fit$df)) + pt(qt(.025,  basin_fit$df), ncp = log(althyp)/slope_se, df = basin_fit$df),
           limit = mean(data$limit),
           first_year = min(data$YEAR),
           n = nrow(data)) %>% 
      mutate(predicted_2018 = exp(intercept + slope * 2018),
             CI = paste0("(", round((exp(confint(basin_fit)[2, 1]) - 1) * 100, 2), ", ", round((exp(confint(basin_fit)[2, 2]) - 1) * 100, 2), ")"),
             ttt = ifelse(is.na(limit), NA, round((log(limit) - intercept) / slope) - 2018),
             ttt = ifelse(ttt < 0, Inf, ttt),
             n.all.lod = sum(data$all.lod))
  }
}

station_row <- function(data, althyp = 1.1){
  station_fit <- lm(log(value) ~ YEAR, data = data)
  tibble(LOC = mocis_name_station(data$LOC[1]),
         order = data$order[1],
         intercept = coef(station_fit)["(Intercept)"],
         slope = coef(station_fit)["YEAR"],
         slope_se = vcov(station_fit)["YEAR", "YEAR"] %>% sqrt(),
         relative_rate = (exp(slope) - 1) * 100,
         p_vs_sat = NA,
         power = (1 - pt(qt(.975,  station_fit$df), ncp = log(althyp)/slope_se, df = station_fit$df)) + pt(qt(.025,  station_fit$df), ncp = log(althyp)/slope_se, df = station_fit$df),
         limit = mean(data$limit),
         first_year = min(data$YEAR),
         n = nrow(data)) %>% 
    mutate(predicted_2018 = exp(intercept + slope * 2018),
           CI = paste0("(", round((exp(confint(station_fit)[2, 1]) - 1) * 100, 2), ", ", round((exp(confint(station_fit)[2, 2]) - 1) * 100, 2), ")"),
           ttt = ifelse(is.na(limit), NA, round((log(limit) - intercept) / slope) - 2018),
           ttt = ifelse(ttt < 0, Inf, ttt),
           n.all.lod = sum(data$all.lod))
}

pretty_basin_kable <- function(table, caption){
  n_col <- ncol(table)
  basin_rows <- which(table$Location == table$basin)
  table %>% select(-basin) %>% 
    mutate_if(is.numeric, round, digits = 3) %>% 
    mutate_all(as.character) %>% 
    mutate_all(~ifelse(is.na(.x), "", .x)) %>% 
    kable(escape = FALSE, digits = 3, caption = caption, align = paste(c("l", rep("r", n_col - 1)), collapse = "")) %>%
    row_spec(basin_rows, bold = TRUE, background = "#E0E0E0") %>% 
    kable_styling()
}

common_intercept <- function(fit.lm){
  theta <- coef(fit.lm)
  d <- length(theta)
  D <- cbind(1, diag(d - 1))
  D[1, 2] <- 0
  intercepts <- t(D %*% theta)
  Sigma <- D %*% vcov(fit.lm) %*% t(D)
  w <- diag(Sigma)
  sum(intercepts/w)/sum(1/w)
}


get_table <- function(data){
  data %>% 
    group_by(basin, var, gen, basin_order) %>% 
    nest() %>% 
    mutate(bas_table = map2(basin, data, basin_table)) %>% 
    select(-data) %>% 
    unnest(cols = bas_table) %>% 
    ungroup()
}

stations <- function(){
  ## Provide a link between the short names used in the data in the labels and
  ## coordinates of the stations.
  list(hav = list("40G7" = list(name = "E Bornholm basin (offsh.)",
                                NKOO = 6181395, EKOO = 1606416),
                  "46H0" = list(name = "N Baltic Proper (offsh.)",
                                NKOO = 6523708, EKOO = 1771651),
                  "49G9" = list(name = "Sea of \U00C5land (offsh.)",
                                NKOO = 6686566, EKOO = 1696248),
                  "51G9" = list(name = "Bothnian Sea (offsh.)",
                                NKOO = 6798326, EKOO = 1698277),
                  ABBE = list(name = "Abbek\U00E5s",
                              NKOO = 6134000, EKOO = 1360700),
                  ANGV = list(name = "\U00C4ngsk\U00E4rsklubb (spring)",
                              NKOO = 6715100, EKOO = 1629400),
                  ANGK = list(name = "\U00C4ngsk\U00E4rsklubb",
                              NKOO = 6715100, EKOO = 1629400),
                  BYXE = list(name = "Byxelkrok",
                              NKOO = 6365800,EKOO = 1571500),
                  FLAD = list(name = "Fladen",
                              NKOO = 6348600, EKOO = 1258800),
                  GAFJ = list(name = "Gaviksfj\U00E4rden",
                              NKOO = 7005100, EKOO = 1642800),
                  HABU = list(name = paste0("W Han\U00F6", "bukten"),
                              NKOO = 6181700, EKOO = 1404600),
                  HAFJ = list(name = "Harufj\U00E4rden",
                              NKOO = 7294000, EKOO = 1825900),
                  HOLM = list(name = paste0("Holm\U00F6","arna"),
                              NKOO = 7073600, EKOO = 1750800),
                  KIFJ = list(name = paste0("Kinnb\U00E4","cksfj\U00E4rden"),
                              NKOO = 7204900, EKOO = 1759200),
                  KULL = list(name = "Kullen",
                              NKOO = 6249400, EKOO = 1288200),
                  LAFJ = list(name = "L\U00E5ngvindsfj\U00E4rden",
                              NKOO = 6852200, EKOO = 1587100),
                  LAGN = list(name = "Lagn\U00F6",
                              NKOO = 6593400, EKOO = 1660100),
                  LAND = list(name = "Landsort",
                              NKOO = 6510000, EKOO = 1627500),
                  RAFJ = list(name = "R\U00E5nefj\U00E4rden",
                              NKOO = 7310900, EKOO = 1802700),
                  UTLA = list(name = "Utl\U00E4ngan",
                              NKOO = 6208830, EKOO = 1501600),
                  UTLV = list(name = "Utl\U00E4ngan (spring)",
                              NKOO = 6208830, EKOO = 1501600),
                  VADO = list(name = paste0("V\U00E4","der\U00F6","arna"),
                              NKOO = 6502000, EKOO = 1218300),
                  SEGO = list(name = "SE Gotland",
                              NKOO = 6294700, EKOO = 1664600),
                  TJRN = list(name = "Tj\U00E4rn\U00F6",
                              NKOO = 6539496, EKOO = 1230900),
                  FJBA = list(name = "Fj\U00E4llbacka",
                              NKOO = 6510100, EKOO = 1236200),
                  KVFJ = list(name = paste0("Kv\U00E4","d\U00F6","fj\U00E4rden"),
                              NKOO = 6434800, EKOO = 1556700),
                  NIDI = list(name = "Nidingen",
                              NKOO = 6368600, EKOO = 1265100),
                  ORFJ = list(name = "\U00D6refj\U00E4rden",
                              NKOO = 7039900, EKOO = 1679300),
                  STKA = list(name = "Stora Karls\U00F6",
                              NKOO = 6352800, EKOO = 1631500)))
}

mocis_name_station <- function(loc){
  map_chr(loc, ~ifelse(.x %in% names(stations()[["hav"]]), 
                       stations()[["hav"]][[.x]][["name"]], paste(.x, "(full name unavailable)")))
}


mocis_get_unit_HTML <- function(var, gen){
  labels <- list(
    fish = list(
      KOND = "Fultons condition factor K",
      FPRC = "Fat % muscle",
      TOTV = "Total weight g",
      TOTL = "Total length cm",
      ALDR = "Age years",
      HBCD = "HBCDD ng/g lw muscle",
      HG = "Hg, ng/g ww muscle",
      PB = "Pb, &mu;g/g dw liver",
      CD = "Cd, &mu;g/g dw liver",
      CU = "Cu, &mu;g/g dw liver",
      ZN = "Zn, &mu;g/g dw liver",
      NI = "Ni, &mu;g/g dw liver",
      CR = "Cr,  &mu;g/g dw liver",
      AG = "Ag,  &mu;g/g dw liver",
      AS = "As &mu;g/g dw liver",
      AL = "Al, &mu;g/g dw liver",
      SE = "Se, &mu;g/g dw liver",
      SN = "Sn, &mu;g/g dw liver",
      DDE = "DDE, &mu;g/g lw muscle",
      DDD = "DDD, &mu;g/g lw muscle",
      DDT = "DDT, &mu;g/g lw muscle",
      PCBSUM = "SigmaPCBs, &mu;g/g lw muscle",
      LINDA = "Lindane, &mu;g/g lw,  muscle",
      AHCH = "alpha-HCH, &mu;g/g lw muscle",
      BHCH = "beta-HCH, &mu;g/g lw muscle",
      HCB = "HCB, &mu;g/g lw muscle",
      CB28 = "PCB-28, &mu;g/g lw muscle",
      CB52 = "PCB-52, &mu;g/g lw muscle",
      CB101 = "PCB-101, &mu;g/g lw muscle",
      CB118 = "PCB-118, &mu;g/g lw muscle",
      CB153 = "PCB-153, &mu;g/g lw muscle",
      CB180 = "PCB-180, &mu;g/g lw muscle",
      BDE47 = "BDE-47, ng/g lw muscle",
      BDE99 = "BDE-99, ng/g lw muscle",
      BDE100 = "BDE-100, ng/g lw muscle",
      BDE153 = "BDE-153, ng/g lw muscle",
      BDE154 = "BDE-154, ng/g lw muscle",
      BDE209 = "BDE-209, ng/g lw muscle",
      BDE28 = "BDE-28, ng/g lw muscle",
      TCDD = "2,3,7,8-TCDD, pg/g lw muscle",
      PECDD = "1,2,3,7,8-PeCDD, pg/g lw muscle",
      HXCDD1 = "1,2,3,4,7,8-HxCDD, pg/g lw muscle",
      HXCDD2 = "1,2,3,6,7,8-HxCDD, pg/g lw muscle",
      HXCDD3 = "1,2,3,7,8,9-HxCDD, pg/g lw muscle",
      HPCDD = "1,2,3,4,6,7,8-HpCDD, pg/g lw muscle",
      OCDD = "1,2,3,4,6,7,8,9-OCDD, pg/g lw muscle",
      TCDF = "2,3,7,8-TCDF, pg/g lw muscle",
      PECDF1 = "1,2,3,7,8-PeCDF, pg/g lw muscle",
      PECDF2 = "2,3,4,7,8-PeCDF, pg/g lw muscle",
      HXCDF1 = "1,2,3,4,7,8-HxCDF, pg/g lw muscle",
      HXCDF2 = "1,2,3,6,7,8-HxCDF, pg/g lw muscle",
      HXCDF3 = "2,3,4,6,7,8-HxCDF, pg/g lw muscle",
      HXCDF4 = "1,2,3,7,8,9-HxCDF, pg/g lw muscle",
      HPCDF1 = "1,2,3,4,6,7,8-HpCDF, pg/g lw muscle",
      HPCDF2 = "1,2,3,4,7,8,9-HpCDF, pg/g lw muscle",
      OCDF = "1,2,3,4,6,7,8,9-OCDF, pg/g lw muscle",
      CB77 = "PCB-77 pg/g lw muscle",
      CB81 = "PCB-81 pg/g lw muscle",
      CB126 = "PCB-126 pg/g lw muscle",
      CB169 = "PCB-169 pg/g lw muscle",
      CB105 = "PCB-105 pg/g lw muscle",
      CB114 = "PCB-114 pg/g lw muscle",
      CB118 = "PCB-118 pg/g lw muscle",
      CB123 = "PCB-123 pg/g lw muscle",
      CB138 = "PCB-138, &mu;g/g lw soft tissue",
      CB156 = "PCB-156 pg/g lw muscle",
      CB157 = "PCB-157 pg/g lw muscle",
      CB167 = "PCB-167 pg/g lw muscle",
      CB189 = "PCB-189 pg/g lw muscle",
      TCDDEQVW = "TEQ WHO -98 (PCDD/F), pg/g ww muscle",
      TCDDEQV = "TEQ WHO -98 (PCDD/F), pg/g lw muscle",
      CBEQV = "TEQ WHO -98 (PCB), pg/g lw muscle",
      TCDDEQ05 = "TEQ WHO -05 (PCDD/F), pg/g lw muscle",
      TCDDEQ05W = "TEQ WHO -05 (PCDD/F), pg/g ww muscle",
      CBEQ05 = "TEQ WHO -05 (PCB), pg/g lw muscle",
      PFHXA = "PFHxA, ng/g ww liver",
      PFHPA = "PFHpA, ng/g ww liver",
      PFOA = "PFOA, ng/g ww liver",
      PFNA = "PFNA, ng/g ww liver",
      PFDA = "PFDA, ng/g ww liver",
      PFUNDA = "PFUnDA, ng/g ww liver",
      PFDODA = "PFDoDA, ng/g ww liver",
      PFTRDA = "PFTrDA, ng/g ww liver",
      PFTEDA = "PFTeTA, ng/g ww liver",
      PFPEDA = "PFPeDA, ng/g ww liver",
      PFBS = "PFBS, ng/g ww liver",
      PFHXS = "PFHxS, ng/g ww liver",
      PFOS = "PFOS, ng/g ww liver",
      PFDS = "PFDS, ng/g ww liver",
      FOSA = "FOSA, ng/g ww liver",
      LFOSA = "lin-FOSA, ng/g ww liver",
      BFOSA = "br-FOSA, ng/g ww liver",
      LPFOS = "lin-PFOS, ng/g  ww liver",
      BPFOS = "br-PFOS, ng/g ww liver",
      LPFDS = "lin-PFDS, ng/g ww liver",
      BPFDS = "br-PFDS, ng/g ww liver",
      D13CUSD = "delta^{13}C,  muscle",
      D15NUSD = "delta^{15}N,  muscle",
      MBT = "Monobytultin, ng/g ww liver",
      DIBT = "Bibutyltin, ng/g ww liver",
      TBT = "Tributyltin (TBT), ng/g ww liver",
      MPT = "Monophenyltin, ng/g ww liver",
      DIPT = "Diphenyltin, ng/g ww liver",
      TPT = "Triphenyltin, ng/g ww liver",
      MOT = "Monooctyltin, ng/g ww liver",
      DIOT = "Dioctyltin, ng/g ww liver)",
      TEQ98lw = "TEQ WHO -98 (PCDD/F+PCB) pg/g lw",
      TEQ98ww = "TEQ WHO -98 (PCDD/F+PCB) pg/g ww",
      TEQ05lw = "TEQ WHO -05 (PCDD/F+PCB) pg/g lw",
      TEQ05ww = "TEQ WHO -05 (PCDD/F+PCB) pg/g ww",      
      PCBsum6 = "Sum of CB-28, 52, 101, 138, 153, 180, &mu;g/g ww",
      PBDEsum5 = "Sum of BDE-47, 99, 100, 153, 154, ng/g ww"),
    bluemussel = list(
      KOND = "Fultons condition factor, K",
      TOTV = "Total weight, g",
      TOTL = "Total length, cm",
      ALDR = "Age, years",
      HG = "Hg, ng/g ww soft tissue",
      PB = "Pb, &mu;g/g ww soft tissue",
      CD = "Cd, &mu;g/g ww soft tissue",
      CU = "Cu, &mu;g/g ww soft tissue",
      ZN = "Zn, &mu;g/g ww soft tissue",
      NI = "Ni, &mu;g/g ww soft tissue",
      CR = "Cr,  &mu;g/g ww soft tissue",
      AG = "Ag,  &mu;g/g ww soft tissue",
      AS = "As, &mu;g/g ww soft tissue",
      AL = "Al, &mu;g/g ww soft tissue",
      SE = "Se, &mu;g/g ww soft tissue",
      SN = "Sn, &mu;g/g ww soft tissue",
      DDE = "DDE, &mu;g/g lw soft tissue",
      DDD = "DDD, &mu;g/g lw soft tissue",
      DDT = "DDT, &mu;g/g lw soft tissue",
      PCBSUM = "SigmaPCBs, &mu;g/g lw soft tissue",
      LINDA = "Lindane, &mu;g/g lw,  soft tissue",
      AHCH = "alpha-HCH, &mu;g/g lw soft tissue",
      BHCH = "beta-HCH, &mu;g/g lw soft tissue",
      HCB = "HCB, &mu;g/g lw soft tissue",
      CB28 = "PCB-28, &mu;g/g lw soft tissue",
      CB52 = "PCB-52, &mu;g/g lw soft tissue",
      CB101 = "PCB-101, &mu;g/g lw soft tissue",
      CB118 = "PCB-118, &mu;g/g lw soft tissue",
      CB138 = "PCB-138, &mu;g/g lw soft tissue",
      CB153 = "PCB-153, &mu;g/g lw soft tissue",
      CB180 = "PCB-180, &mu;g/g lw soft tissue",
      BDE47 = "BDE-47, ng/g lw soft tissue",
      BDE99 = "BDE-99, ng/g lw soft tissue",
      BDE100 = "BDE-100, ng/g lw soft tissue",
      BDE153 = "BDE-153, ng/g lw soft tissue",
      BDE154 = "BDE-154, ng/g lw soft tissue",
      BDE209 = "BDE-209, ng/g lw soft tissue",
      BDE28 = "BDE-28, ng/g lw soft tissue",
      NAP = "Naphtalene, ng/g dw soft tissue",
      ACNE = "Acenapthene, ng/g dw soft tissue",
      FLE = "Fluorene, ng/g dw soft tissue",
      PA = "Phenantrene, ng/g dw soft tissue",
      ANT = "Anthracene, ng/g dw soft tissue",
      FLU = "Fluoranthene, ng/g dw soft tissue",
      PYR = "Pyrene, ng/g dw soft tissue",
      BAA = "Benso(a)anthracene, ng/g dw soft tissue",
      CHR = "Chrysene, ng/g dw soft tissue",
      BBF = "Benso(b)fluoranthene, ng/g dw soft tissue",
      BKF = "Benso(k)fluoranthene, ng/g dw soft tissue",
      BAP = "Benso(a)pyrene, ng/g dw soft tissue",
      DBAHA = "Dibenso(a,h)anthracene, ng/g dw soft tissue",
      BGHIP = "Benso(g,h,i)perylene, ng/g dw soft tissue",
      ICDP = "Indeno(1,2,3-cd)pyrene, ng/g dw soft tissue",
      SUMPAH = "SigmaPAHs, ng/g dw soft tissue",
      HBCD = "HBCDD, ng/g dw soft tissue",
      D13CUSD = "delta^{13}C,  soft tissue",
      D15NUSD = "delta^{15}N,  soft tissue)",
      TEQ98lw = "TEQ WHO -98 (PCDD/F+PCB) pg/g lw",
      TEQ98ww = "TEQ WHO -98 (PCDD/F+PCB) pg/g ww",
      TEQ05lw = "TEQ WHO -05 (PCDD/F+PCB) pg/g lw",
      TEQ05ww = "TEQ WHO -05 (PCDD/F+PCB) pg/g ww",   
      PCBsum6 = "Sum of CB-28, 52, 101, 138, 153, 180, &mu;g/g ww",
      PBDEsum5 = "Sum of BDE-47, 99, 100, 153, 154, ng/g ww"),
    birds = list(
      HG = "Hg, ng/g ww egg",
      PB = "Pb, &mu;g/g dw egg",
      CD = "Cd, &mu;g/g dw egg",
      CU = "Cu, &mu;g/g dw egg",
      ZN = "Zn, &mu;g/g dw egg",
      NI = "Ni, &mu;g/g dw egg",
      CR = "Cr,  &mu;g/g dw egg",
      AG = "Ag,  &mu;g/g dw egg",
      AS = "As, &mu;g/g dw egg",
      AL = "Al, &mu;g/g dw egg",
      SE = "Se, &mu;g/g dw egg",
      SN = "Sn, &mu;g/g dw egg",
      DDE = "DDE, &mu;g/g lw egg",
      DDD = "DDD, &mu;g/g lw egg",
      DDT = "DDT, &mu;g/g lw egg",
      PCBSUM = "SigmaPCBs, &mu;g/g lw egg",
      LINDA = "Lindane, &mu;g/g lw,  egg",
      AHCH = "alpha-HCH, &mu;g/g lw egg",
      BHCH = "beta-HCH, &mu;g/g lw egg",
      HCB = "HCB, &mu;g/g lw egg",
      CB28 = "PCB-28, &mu;g/g lw egg",
      CB52 = "PCB-52, &mu;g/g lw egg",
      CB101 = "PCB-101, &mu;g/g lw egg",
      CB118 = "PCB-118, &mu;g/g lw egg",
      CB138 = "PCB-138, &mu;g/g lw soft tissue",
      CB153 = "PCB-153, &mu;g/g lw egg",
      CB180 = "PCB-180, &mu;g/g lw egg",
      BDE47 = "BDE-47, ng/g lw egg",
      BDE99 = "BDE-99, ng/g lw egg",
      BDE100 = "BDE-100, ng/g lw egg",
      BDE153 = "BDE-153, ng/g lw egg",
      BDE154 = "BDE-154, ng/g lw egg",
      BDE209 = "BDE-209, ng/g lw egg",
      BDE28 = "BDE-28, ng/g lw egg",
      TCDD = "2,3,7,8-TCDD, pg/g lw egg",
      PECDD = "1,2,3,7,8-PeCDD, pg/g lw egg",
      HXCDD1 = "1,2,3,4,7,8-HxCDD, pg/g lw egg",
      HXCDD2 = "1,2,3,6,7,8-HxCDD, pg/g lw egg",
      HXCDD3 = "1,2,3,7,8,9-HxCDD, pg/g lw egg",
      HPCDD = "1,2,3,4,6,7,8-HpCDD, pg/g lw egg",
      OCDD = "1,2,3,4,6,7,8,9-OCDD, pg/g lw egg",
      TCDF = "2,3,7,8-TCDF, pg/g lw egg",
      PECDF1 = "1,2,3,7,8-PeCDF, pg/g lw egg",
      PECDF2 = "2,3,4,7,8-PeCDF, pg/g lw egg",
      HXCDF1 = "1,2,3,4,7,8-HxCDF, pg/g lw egg",
      HXCDF2 = "1,2,3,6,7,8-HxCDF, pg/g lw egg",
      HXCDF3 = "2,3,4,6,7,8-HxCDF, pg/g lw egg",
      HXCDF4 = "1,2,3,7,8,9-HxCDF, pg/g lw egg",
      HPCDF1 = "1,2,3,4,6,7,8-HpCDF, pg/g lw egg",
      HPCDF2 = "1,2,3,4,7,8,9-HpCDF, pg/g lw egg",
      OCDF = "1,2,3,4,6,7,8,9-OCDF, pg/g lw egg",
      CB77 = "PCB-77 pg/g lw egg",
      CB81 = "PCB-81 pg/g lw egg",
      CB126 = "PCB-126 pg/g lw egg",
      CB169 = "PCB-169 pg/g lw egg",
      CB105 = "PCB-105 pg/g lw egg",
      CB114 = "PCB-114 pg/g lw egg",
      CB118 = "PCB-118 pg/g lw egg",
      CB123 = "PCB-123 pg/g lw egg",
      CB156 = "PCB-156 pg/g lw egg",
      CB157 = "PCB-157 pg/g lw egg",
      CB167 = "PCB-167 pg/g lw egg",
      CB189 = "PCB-189 pg/g lw egg",
      TCDDEQVW = "TEQ WHO -98 (PCDD/F), pg/g ww egg",
      TCDDEQV = "TEQ WHO -98 (PCDD/F), pg/g lw egg",
      CBEQV = "TEQ WHO -98 (PCB), pg/g lw egg",
      TCDDEQ05 = "TEQ WHO -05 (PCDD/F), pg/g lw egg",
      TCDDEQ05W = "TEQ WHO -05 (PCDD/F), pg/g ww egg",
      CBEQ05 = "TEQ WHO -05 (PCB), pg/g lw egg",
      PFHXA = "PFHxA, ng/g ww egg",
      PFHPA = "PFHpA, ng/g ww egg",
      PFOA = "PFOA, ng/g ww egg",
      PFNA = "PFNA, ng/g ww egg",
      PFDA = "PFDA, ng/g ww egg",
      PFUNDA = "PFUnDA, ng/g ww egg",
      PFDODA = "PFDoDA, ng/g ww egg",
      PFTRDA = "PFTrDA, ng/g ww egg",
      PFTEDA = "PFTeTA, ng/g ww egg",
      PFPEDA = "PFPeDA, ng/g ww egg",
      PFBS = "PFBS, ng/g ww egg",
      PFHXS = "PFHxS, ng/g ww egg",
      PFOS = "PFOS, ng/g ww egg",
      PFDS = "PFDS, ng/g ww egg",
      FOSA = "FOSA, ng/g ww egg",
      LFOSA = "lin-FOSA, ng/g ww egg",
      BFOSA = "br-FOSA, ng/g ww egg",
      LPFOS = "lin-PFOS, ng/g  ww egg",
      BPFOS = "br-PFOS, ng/g ww egg",
      LPFDS = "lin-PFDS, ng/g ww egg",
      BPFDS = "br-PFDS, ng/g ww egg",
      HBCD = "HBCDD, ng/g ww egg",
      D13CUSD = "delta^{13}C,  egg",
      D15NUSD = "delta^{15}N,  egg",
      TEQ98lw = "TEQ WHO -98 (PCDD/F+PCB) pg/g lw",
      TEQ98ww = "TEQ WHO -98 (PCDD/F+PCB) pg/g ww",
      TEQ05lw = "TEQ WHO -05 (PCDD/F+PCB) pg/g lw",
      TEQ05ww = "TEQ WHO -05 (PCDD/F+PCB) pg/g ww",   
      PCBsum6 = "Sum of CB-28, 52, 101, 138, 153, 180, &mu;g/g ww",
      PBDEsum5 = "Sum of BDE-47, 99, 100, 153, 154, ng/g ww"))
  label <- ""
  if (gen %in% c("Herring", "Perch", "Cod", "Eelpout", "Pike", "Arctic char")) {
    label <- getElement(getElement(labels, "fish"), var)
  }
  if (gen == "Blue mussel") {
    label <- getElement(getElement(labels, "bluemussel"), var)
  }
  if (gen %in% c("Guillemot", "Eurasian Oystercatcher", "Common tern")) {
    label <- getElement(getElement(labels, "birds"), var)
  }
  ifelse(is.null(label), "", label)
}