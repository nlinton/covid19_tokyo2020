libs = c('dplyr','magrittr','tidyr','purrr','data.table','readxl','ggplot2','cowplot')
for(x in libs) { suppressMessages(library(x, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE)) }

source("bp_func.R")

# Updated for R version 4.1.1




### DATA

load("bp_data.rda")
# OWID data obtained from https://covid.ourworldindata.org/data/
# Japan and Tokyo case data obtained from https://covid19.mhlw.go.jp/
# Tokyo 2020 case and testing data available at https://github.com/nlinton/covid19_tokyo2020
# Tokyo vaccination data obtained from https://cio.go.jp/c19vaccine_opendata on 2021-07-15




### DATES

# Olympics dates
gen_length = 5 # based on serial interval of wild-type covid-19
olympics_start = as.Date("2021-07-21")
olympics_opening = as.Date("2021-07-23")
olympics_end = as.Date("2021-08-08")

length_of_olympics = as.numeric(olympics_end - olympics_start + 1) ; sprintf("The Olympics last %i days", length_of_olympics)
gen_in_olympics = ceiling((length_of_olympics)/gen_length)
G <- gen_in_olympics + 2 # one before, one after

analysis_start = olympics_start - gen_length
analysis_end = analysis_start + (G*gen_length) - 1

data_cutoff = olympics_start - gen_length - 1 ; data_cutoff
today <- Sys.Date()





### PARAMETERS AND VARIABLES

# Sensitivity parameters
grp_names_3 <- c("Games accredited", "Tokyo population", "Domestic spectators")
grp_names_2 <- c("Games accredited", "Tokyo population")
tokyo_z0 = c(1200,1500,1800)*gen_length
Rb = c(0.7, 0.9, 1.2)
k = c(0.2, 0.6)
theta = c(0.95, 0.99)
G = G
nsim=10000
loop_grid <- expand.grid(Rb, k, tokyo_z0/gen_length, theta) |> rename(Rb=1, k=2, tokyo_inc=3, theta=4)

# Vaccine parameters
arr_vax = 0.8 # % of foreign arrivals vaccinated
antigen_sens = 0.5 # sensitivity of antigen tests done at airport (on the low end)
ve_voc_delta_1dose = 0.3
ve_voc_delta_2dose = 0.8
world_ve = 0.6 # Considers some vaccines are less effective than others






### DETERMINE CASE PREVALENCES

prev_length <- 7 # window for determining case prevalence

# World
world_pop <- (owid |> filter(location=="World" & date==data_cutoff))[,"population"] ; sprintf("World population: %.0f", world_pop)
world_tail <- tail(owid |> filter(location=="World" & date<data_cutoff), prev_length)[,"new_cases"]
world_recovered <- 173175848 # 7/15 recovered, https://www.worldometers.info/coronavirus/coronavirus-cases/
world_anydose <- (owid |> filter(location=="World" & date==data_cutoff))[,"people_vaccinated"]
world_prev <- sum(world_tail)/(world_pop - world_anydose*world_ve - world_recovered)

# Japan
japan_pop <- (owid |> filter(location=="Japan" & date==data_cutoff))[,"population"] ; sprintf("Japan population: %.0f", japan_pop)
japan_tail <- tail(japan_cases |> filter(date<data_cutoff), prev_length)[,"npos"]
japan_recovered = 791751 # 7/15 recovered, https://www.mhlw.go.jp/stf/newpage_19871.html
japan_vax <- (owid |> filter(location=="Japan" & date==data_cutoff))[,c("date","location","people_vaccinated","people_fully_vaccinated")]
japan_2dose <- japan_vax[,"people_fully_vaccinated"] ; sprintf("Japan 1-dose: %.0f", japan_2dose)
japan_1dose <- japan_vax[, "people_vaccinated"] - japan_2dose ; sprintf("Japan 2-dose: %.0f", japan_1dose)
japan_frac_2dose <- japan_2dose/japan_pop
japan_frac_1dose <- japan_1dose/japan_pop
japan_prev <- sum(japan_tail)/(japan_pop - japan_1dose*ve_voc_delta_1dose - japan_2dose*ve_voc_delta_2dose - japan_recovered)


# Tokyo
tokyo_pop = 13960000 # Approximate Tokyo population
tokyo_vax <- jpn_vax_vrs |> filter(prefecture=="13") |> group_by(status) |> summarize(tot=sum(count))
tokyo_2dose <- as.numeric(tokyo_vax[2,2])
tokyo_1dose <- as.numeric(tokyo_vax[1,2])
tokyo_frac_1dose = tokyo_1dose/tokyo_pop; sprintf("Fraction of Tokyo population vaccinated with 1 dose %.3f", tokyo_frac_1dose)
tokyo_frac_2dose = tokyo_2dose/tokyo_pop; sprintf("Fraction of Tokyo population vaccinated with 2 doses %.3f", tokyo_frac_2dose)
tokyo_recovered = 179473 # 7/20, https://stopcovid19.metro.tokyo.lg.jp/
tokyo_prev = tokyo_z0/(tokyo_pop - tokyo_2dose*ve_voc_delta_2dose - tokyo_1dose*ve_voc_delta_1dose - tokyo_recovered)







### Determine number of index cases by group


### Olympics

# Athletes and foreign participants
games_foreign_athletes <- 11000
games_foreign_other <- 31000
games_foreign <- games_foreign_athletes + games_foreign_other


games_foreign_cases = (games_foreign*(1-antigen_sens)*(1-arr_vax*ve_voc_delta_2dose)*world_prev)

# Domestic participants
reduction_in_volunteers = 30000
frac_inv = 0.8
games_volunteers = 54000
games_volunteers_red = games_volunteers - reduction_in_volunteers
games_domestic = 190000
games_domestic_red = games_domestic - reduction_in_volunteers
games_all = games_domestic + games_foreign
games_all_red = games_all - games_volunteers_red
games_domestic_cases = (games_domestic*(1-tokyo_frac_2dose*ve_voc_delta_2dose)*(1-tokyo_frac_1dose*ve_voc_delta_1dose)*tokyo_prev)
games_domestic_cases_red = (games_domestic_red*(1-tokyo_frac_2dose*ve_voc_delta_2dose)*(1-tokyo_frac_1dose*ve_voc_delta_1dose)*tokyo_prev)


# All games participants cases
games_cases <- ceiling(games_foreign_cases + games_domestic_cases)
games_cases_red <- ceiling(games_foreign_cases + games_domestic_cases_red)
games_cases
games_cases_red
games_foreign_cases + games_domestic_cases


# All games participants

# Distribute between generations by percentage
distribute_gp_z0 <- c(0.4,0.4,0.1,0.1,0,0)


# For scenario with spectators
npart_gen_o <- cbind(expand.grid(gen=1:G-1, type=grp_names_3[1], tokyo_inc=tokyo_z0/gen_length, Rb=Rb, k=k, theta=theta),
                     n=rep(ceiling(c((games_all-games_volunteers)*frac_inv,
                                     rep((games_all-games_volunteers*.5)*frac_inv, G-1))), length(tokyo_z0)),  # No volunteers in gen 0,
                     z0=rep(c(ceiling(sapply(1:length(tokyo_z0), \(j) games_cases[j] * distribute_gp_z0))),length(tokyo_z0)))
npart_gen_o

# For scenario without spectators
npart_gen_o_red <- cbind(expand.grid(gen=1:G-1, type=grp_names_3[1], tokyo_inc=tokyo_z0/gen_length, Rb=Rb, k=k, theta=theta),
                         n=rep(ceiling(c((games_all_red-games_volunteers_red)*frac_inv,
                                         rep((games_all_red-games_volunteers_red*.5)*frac_inv, G-1))), length(tokyo_z0)),  # reduce volunteers
                         z0=rep(c(ceiling(sapply(1:length(tokyo_z0), \(j) games_cases_red[j] * distribute_gp_z0))),length(tokyo_z0)))
npart_gen_o_red

# Tokyo
ntokyo_gen_o <- cbind(expand.grid(gen=1:G-1, type=grp_names_3[2], tokyo_inc=tokyo_z0/gen_length, Rb=Rb, k=k, theta=theta),
                      n=tokyo_pop,
                      z0=as.numeric(sapply(seq_along(grp_names_3), \(j) c(tokyo_z0[j], rep(0,G-1)))))
rbind(ntokyo_gen_o |> head(G), ntokyo_gen_o |> tail(G))


# Total spectators by generation
daily_spec_tokyo_o <- venues |> filter(admin1=="東京" & date<olympics_end+1) |> group_by(date) |> summarize(nspec_daily=sum(daily_occ, na.rm=T))
spec_gen_o <- cbind(expand.grid(gen=1:G-1, type=grp_names_3[3], tokyo_inc=tokyo_z0/gen_length, Rb=Rb, k=k, theta=theta)) |>
  left_join(rbind(c(0,0),daily_spec_tokyo_o |> mutate(gen=rep(1:(G-1), each=gen_length)[1:nrow(daily_spec_tokyo_o)]) |> group_by(gen) |> summarize(n=sum(nspec_daily))), by="gen") |>
  mutate(z0=round(n*(1-japan_frac_2dose*ve_voc_delta_2dose)*(1-japan_frac_1dose*ve_voc_delta_1dose)*japan_prev,0)) |>
  replace_na(list(n=0,z0=0))





#############################################



### MODEL OUTPUT



### 3 groups scenario ###


# Next-generation matrix for three groups analysis
D <- length(grp_names_3) # Set number of groups
df3 <- rbind(npart_gen_o, ntokyo_gen_o, spec_gen_o)

df_for_ngm <- df3 |> dplyr::select(gen, type, n) |> unique()

ngm3 <- lapply(seq_along(Rb), \(b) {
  lapply(seq_along(theta), \(t) {
    create_ngm(df_for_ngm, G, D, Rb[b], theta[t], grp_names_3) }) })

names(ngm3) <- paste0("Rb_", Rb) ; names(ngm3)
ngm3 <- lapply(ngm3, \(x) { names(x) <- paste0("theta_", theta) ; x } )



# Initial cases
z0_ngm3 <- lapply(1:3, \(z) matrix(df3 |> dplyr::select(gen, type, tokyo_inc, z0) |> unique() |>
                                     filter(tokyo_inc==tokyo_z0[z]/gen_length) |> dplyr::select(z0) |> pull(), nrow=G))
names(z0_ngm3) <- paste0("tokyo_", tokyo_z0/gen_length) # z0 only varies by tokyo_inc



# 3 groups with total output per group per generation
res3_3d <- lapply(1:nrow(loop_grid), \(l) {
  df <- loop_grid[l,]
  Rb_ <- paste0("Rb_", df$Rb)
  k_ <- df$k
  theta_ <- paste0("theta_", df$theta)
  tokyo_inc_ <- paste0("tokyo_", df$tokyo_inc)
  z0_ <- z0_ngm3[[tokyo_inc_]]
  ngm_ <- ngm3[[Rb_]][[theta_]]
  array(unlist(bp_sim(nsim, G, D, z0_, ngm_, k_)), c(G, D, nsim))
})



# 3 groups with total output ij cell per generation
res3_4d <- lapply(1:nrow(loop_grid), \(l) {
  df <- loop_grid[l,]
  Rb_ <- paste0("Rb_", df$Rb)
  k_ <- df$k
  theta_ <- paste0("theta_", df$theta)
  tokyo_inc_ <- paste0("tokyo_", df$tokyo_inc)
  z0_ <- z0_ngm3[[tokyo_inc_]]
  ngm_ <- ngm3[[Rb_]][[theta_]]
  array(unlist(bp_sim_ij(nsim, G, D, z0_, ngm_, k_)),
        c(D, D, G, nsim),
        dimnames=list(i=grp_names_3, j=grp_names_3, gen=(1:G)-1, sim=1:nsim))
})



### 2 groups scenario ###


# Next-generation matrix for three groups analysis
D <- length(grp_names_2)
df2 <- rbind(npart_gen_o, ntokyo_gen_o)

df_for_ngm <- df2 |> dplyr::select(gen, type, n) |> unique()

ngm2 <- lapply(seq_along(Rb), \(b) {
  lapply(seq_along(theta), \(t) {
    create_ngm(df_for_ngm, G, D, Rb[b], theta[t], grp_names_2) }) })

names(ngm2) <- paste0("Rb_", Rb) ; names(ngm2)
ngm2 <- lapply(ngm2, \(x) { names(x) <- paste0("theta_", theta) ; x } )



# Initial cases
z0_ngm2 <- lapply(1:3, \(z) matrix(df2 |> dplyr::select(gen, type, tokyo_inc, z0) |> unique() |>
                                     filter(tokyo_inc==tokyo_z0[z]/5) |> dplyr::select(z0) |> pull(), nrow=G))
names(z0_ngm2) <- paste0("tokyo_", tokyo_z0/gen_length) # z0 only varies by tokyo_inc




# 2 groups with total output per group per generation
res2_3d <- lapply(1:nrow(loop_grid), \(l) {
  df <- loop_grid[l,]
  Rb_ <- paste0("Rb_", df$Rb)
  k_ <- df$k
  theta_ <- paste0("theta_", df$theta)
  tokyo_inc_ <- paste0("tokyo_", df$tokyo_inc)
  z0_ <- z0_ngm2[[tokyo_inc_]]
  ngm_ <- ngm2[[Rb_]][[theta_]]
  array(unlist(bp_sim(nsim, G, D, z0_, ngm_, k_)), c(G, D, nsim))
})




# 2 groups with total output ij cell per generation
res2_4d <- lapply(1:nrow(loop_grid), \(l) {
  df <- loop_grid[l,]
  Rb_ <- paste0("Rb_", df$Rb)
  k_ <- df$k
  theta_ <- paste0("theta_", df$theta)
  tokyo_inc_ <- paste0("tokyo_", df$tokyo_inc)
  z0_ <- z0_ngm2[[tokyo_inc_]]
  ngm_ <- ngm2[[Rb_]][[theta_]]
  array(unlist(bp_sim_ij(nsim, G, D, z0_, ngm_, k_)),
        c(D, D, G, nsim),
        dimnames=list(i=grp_names_2, j=grp_names_2, gen=1:G-1, sim=1:nsim))
})
