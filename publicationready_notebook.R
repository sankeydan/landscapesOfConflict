# clear Work Space
rm ( list = ls( ))

# libraries
{
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(multcomp)
  library(DHARMa)
  library(glmmTMB)
  library(car)
  library(viridis)
}

# load data
load ( file.path ( getwd() ,  "dat_den2023_08_14.rda"))

# descriptive stats
{
  table ( dat_den$igi_TF)
  nrow ( dat_den ) - length ( which ( is.na ( dat_den$riv.max.fami )))
  table ( dat_den$i.l.i_focal)
  mean(table ( dat_den$i.l.i_focal))
  dat_den$coll.hist.both = NA
  dat_den$coll.hist.both[ dat_den$foc.prop.coll.vs.hist == 0] = 0
  dat_den$coll.hist.both[ dat_den$foc.prop.coll.vs.hist == 1] = 1
  dat_den$coll.hist.both[ dat_den$foc.prop.coll.vs.hist < 1 &
                            dat_den$foc.prop.coll.vs.hist > 0] = 0.5
  table(dat_den$coll.hist.both)
  par(mfrow= c(2,1))
  summary ( dat_den$terr.size[dat_den$bby_TF==T] )
  summary ( dat_den$terr.size[dat_den$bby_TF==F] )
  summary ( dat_den$terr.size[dat_den$foc.oes==T])
  summary ( dat_den$terr.size[dat_den$foc.oes==F])
  summary ( dat_den$terr.size[dat_den$igi_TF==T])
  summary ( dat_den$terr.size[dat_den$igi_TF==F])
  summary ( dat_den$terr.size)
}

##### MODEL 1 #####
runmodel1 = T
if(runmodel1){
  # models



  {
    IGIYN_plus <- glmmTMB( foc.fami~
                             bby_TF +
                             igi_TF
                           + scale ( rain_md)
                           +  foc.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           , data= dat_den
                           ,
                           beta_family(link = "logit")
    )
    print("1/8")

    IGIYN_times <- glmmTMB( foc.fami~
                              bby_TF *
                              igi_TF
                            + scale ( rain_md)
                            +  foc.prop.coll.vs.hist
                            + comp
                            + (1| i.l.i_focal)
                            + ar1( datenF-1|acor)
                            , data= dat_den
                            ,
                            beta_family(link = "logit")
    )
    print("2/8")

    IGIFREQ_plus <- glmmTMB( foc.fami~
                               bby_TF +
                               igi_freq
                             + scale ( rain_md)
                             +  foc.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("3/8")

    IGIFREQ_times <- glmmTMB( foc.fami ~
                                bby_TF *
                                igi_freq
                              + scale ( rain_md)
                              + foc.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("4/8")

    nullmodel <- glmmTMB( foc.fami~
                            scale ( rain_md)
                          + foc.prop.coll.vs.hist
                          + comp
                          + (1| i.l.i_focal)
                          + ar1( datenF-1|acor)
                          , data= dat_den
                          ,
                          beta_family(link = "logit")
    )
    print("5/8")

    IGIFREQ_only <- glmmTMB( foc.fami ~
                             #   bby_TF *
                                igi_freq
                              + scale ( rain_md)
                              + foc.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("6/8")

    IGITF_only <- glmmTMB( foc.fami ~
                               #   bby_TF *
                               igi_TF
                             + scale ( rain_md)
                             + foc.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("7/8")


    BABY_only <- glmmTMB( foc.fami ~
                              bby_TF
                             #igi_TF
                           + scale ( rain_md)
                           + foc.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           , data= dat_den
                           ,
                           beta_family(link = "logit")
    )
    print("8/8")
  }

  # AIC
  AICtab = AIC(nullmodel,
               IGIYN_times,
               IGIFREQ_times,
               IGIYN_plus,
               IGIFREQ_plus,
               IGIFREQ_only,
               IGITF_only,
               BABY_only)
  AICtab = AICtab[order(AICtab$AIC),]
  AIC= cbind("Focal familiarity BBY", dimnames(AICtab)[[1]],AICtab,c(AICtab$AIC-AICtab$AIC[1])) # combined is the best AIC
  names(AIC) = c( "Model", "Covariates" , "df" , "AIC" , "deltaAIC")
  AIC


  # diagnostics
  mod=IGIYN_plus
  simulateResiduals(mod,plot= T)
  acf(residuals(mod))


  # statistics
  sum = summary(mod)
  stats = sum$coefficients$cond
  stats = apply (stats , 2 , function (x)  round (x,4))
  stats = cbind ( "foc.fami BBY", stats)
  stats

  # save
  MOD1 = list ( AIC = AIC, stats = stats)
  save(MOD1, file= file.path(getwd() ,  "mod1.rda"))
} else {
  load ( file.path(getwd() ,  "mod1.rda"))
}

# stats
MOD1$AIC
MOD1$stats

# plot
par(mfrow= c(1,2))
boxplot ( dat_den$foc.fami ~ dat_den$bby_TF)
boxplot ( dat_den$foc.fami ~ dat_den$igi_TF)



###### MODEL 2 ######
runmodel2 = T
if(runmodel2){
  # models

  {
    IGITF_plus <- glmmTMB( riv.max.fami^0.2
                           ~
                             bby_TF+
                             #foc.oes+
                             igi_TF
                           + scale ( rain_md)
                           +  riv.max.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           ,data= dat_den
                           ,
                           beta_family(link = "logit")

    )
    print("1/8")

    IGITF_times <- glmmTMB( riv.max.fami^0.2
                            ~
                              bby_TF*
                              igi_TF
                            + scale ( rain_md)
                            +  riv.max.prop.coll.vs.hist
                            + comp
                            + (1| i.l.i_focal)
                            + ar1( datenF-1|acor)
                            , data= dat_den
                            ,
                            beta_family(link = "logit")
    )
    print("2/8")

    IGIFREQ_plus <- glmmTMB( riv.max.fami^0.2
                             ~
                               #foc.oes+
                               bby_TF+
                               riv.max.igi_freq
                             + scale ( rain_md)
                             +  riv.max.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("3/8")

    IGIFREQ_times <- glmmTMB( riv.max.fami^0.2
                              ~
                                bby_TF*
                                riv.max.igi_freq
                              + scale ( rain_md)
                              +  riv.max.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("4/8")

    nullmodel <- glmmTMB( riv.max.fami^0.2
                          ~
                            scale ( rain_md)
                          +  riv.max.prop.coll.vs.hist
                          + comp
                          + (1| i.l.i_focal)
                          + ar1( datenF-1|acor)
                          , data= dat_den
                          ,
                          beta_family(link = "logit")
    )
    print("5/8")

    IGIFREQ_only <- glmmTMB( riv.max.fami^0.2
                              ~
                                #bby_TF*
                                riv.max.igi_freq
                              + scale ( rain_md)
                              +  riv.max.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("6/8")

    IGITF_only <- glmmTMB( riv.max.fami^0.2
                             ~
                               #bby_TF*
                               igi_TF
                             + scale ( rain_md)
                             +  riv.max.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("7/8")

    BABY_only <- glmmTMB( riv.max.fami^0.2
                           ~
                             bby_TF
                             #igi_TF
                           + scale ( rain_md)
                           +  riv.max.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           , data= dat_den
                           ,
                           beta_family(link = "logit")
    )
    print("8/8")
  }

  table(dat_den$bby_TF,dat_den$foc.oes)

  # AIC
  AICtab = AIC(nullmodel,
               IGITF_times,
               IGIFREQ_times,
               IGITF_plus,
               IGIFREQ_plus,
                IGIFREQ_only,
               IGITF_only,
               BABY_only)
  AICtab
  AICtab = AICtab[order(AICtab$AIC),]
  AIC= cbind("Rival familiarity BBY", dimnames(AICtab)[[1]],AICtab,c(AICtab$AIC-AICtab$AIC[1])) # combined is the best AIC
  names(AIC) = c( "Model", "Covariates" , "df" , "AIC" , "deltaAIC")
  AIC

  #diagnostics
  mod = IGIFREQ_plus
  simulateResiduals(mod,plot=T)# not great
  hist(residuals  (mod)) # doesn't look too bad
  qqnorm(residuals(mod))# so try normal qq plot
  qqline(residuals(mod)) # looks fine

  # statistics
  sum = summary(mod)
  stats =sum$coefficients$cond
  stats=apply (stats , 2 , function (x)  round (x,4))
  stats = cbind ( "riv.fami BBY", stats)
  stats


  # save
  MOD2 = list ( AIC = AIC, stats = stats)
  save(MOD2, file= file.path(getwd() ,  "mod2.rda"))
} else {
  load ( file.path(getwd() ,  "mod2.rda"))
}

# stats
MOD2$AIC
MOD2$stats

# plot
cols = viridis(10)
g1 = ggplot(dat_den, aes(y = riv.max.fami, x = riv.max.igi_freq * 90, color = bbyst_Y_N)) +
  geom_point() +
  scale_color_manual(values = c("#482878FF", "#1F9E89FF")) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  coord_cartesian(ylim = c(0,0.954)) +
  labs(color = "Babysitting (Y/N)", y = "Outgroup Familiarity Index (1-CDF)", x = "Frequency of Interactions with Outgroup (in 90 days)") +
  theme(
    axis.title.x = element_text(size = 14),  # adjust x-axis label
    axis.title.y = element_text(size = 14),  # adjust y-axis label
    axis.text.x = element_text(size = 12),   # adjust x-axis ticks text
    axis.text.y = element_text(size = 12),   # adjust y-axis ticks text
    legend.title = element_text(size = 14),  # adjust legend title
    legend.text = element_text(size = 12)    # adjust legend text
  )
g2 = ggplot(dat_den, aes(y = riv.max.fami.subt, x = riv.max.igi_freq * 90, color = bbyst_Y_N)) +
  geom_point() +
  scale_color_manual(values = c("#482878FF", "#1F9E89FF")) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(color = "Babysitting (Y/N)", y = "Familiarity Difference (Focal - Rival)", x = "Frequency of Interactions with Outgroup (in 90 days)") +
  theme(

    axis.title.x = element_text(size = 14),  # adjust x-axis label
    axis.title.y = element_text(size = 14),  # adjust y-axis label
    axis.text.x = element_text(size = 12),   # adjust x-axis ticks text
    axis.text.y = element_text(size = 12),   # adjust y-axis ticks text
    legend.title = element_text(size = 14),  # adjust legend title
    legend.text = element_text(size = 12)    # adjust legend text
  )
grid.arrange( g1, g2,nrow=  1)

#### MODEL 3 #### OESTRUS
runmodel3 = T
if(runmodel3){
  # models
  #dat_den = dat_den[ - which (  dat_den$bby_TF == T & dat_den$foc.oes == T),]
  {
    IGIYN_plus <- glmmTMB( foc.fami~
                             foc.oes +
                             igi_TF
                           + scale ( rain_md)
                           +  foc.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           , data= dat_den
                           ,
                           beta_family(link = "logit")
    )
    print("1/8")

    IGIYN_times <- glmmTMB( foc.fami~
                              foc.oes *
                              igi_TF
                            + scale ( rain_md)
                            +  foc.prop.coll.vs.hist
                            + comp
                            + (1| i.l.i_focal)
                            + ar1( datenF-1|acor)
                            , data= dat_den
                            ,
                            beta_family(link = "logit")
    )
    print("2/8")

    IGIFREQ_plus <- glmmTMB( foc.fami~
                               foc.oes +
                               igi_freq
                             + scale ( rain_md)
                             +  foc.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("3/8")

    IGIFREQ_times <- glmmTMB( foc.fami ~
                                foc.oes *
                                igi_freq
                              + scale ( rain_md)
                              + foc.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("4/8")

    nullmodel <- glmmTMB( foc.fami~
                            scale ( rain_md)
                          + foc.prop.coll.vs.hist
                          + comp
                          + (1| i.l.i_focal)
                          + ar1( datenF-1|acor)
                          , data= dat_den
                          ,
                          beta_family(link = "logit")
    )
    print("5/8")

    IGIFREQ_only <- glmmTMB( foc.fami ~
                            #    foc.oes *
                                igi_freq
                              + scale ( rain_md)
                              + foc.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("6/8")

    IGITF_only <- glmmTMB( foc.fami ~
                               #    foc.oes *
                               igi_TF
                             + scale ( rain_md)
                             + foc.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("7/8")

    OES_only <- glmmTMB( foc.fami ~
                             foc.oes

                           + scale ( rain_md)
                           + foc.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           , data= dat_den
                           ,
                           beta_family(link = "logit")
    )
    print("8/8")
  }

  # # bonus model #  #
  IGIYN_plusOESandBBY <- glmmTMB( foc.fami~
                           foc.oes +
                           igi_TF +
                           bby_TF
                         + scale ( rain_md)
                         +  foc.prop.coll.vs.hist
                         + comp
                         + (1| i.l.i_focal)
                         + ar1( datenF-1|acor)
                         , data= dat_den
                         ,
                         beta_family(link = "logit")
  )

  summary(IGIYN_plusOESandBBY)
  # AIC
  AICtab = AIC(nullmodel,
               IGIYN_times,
               IGIFREQ_times,
               IGIYN_plus,
               IGIFREQ_plus,
               IGIFREQ_only,
               IGITF_only,
               OES_only)
  AICtab = AICtab[order(AICtab$AIC),]
  AIC= cbind("Focal familiarity OES", dimnames(AICtab)[[1]],AICtab,c(AICtab$AIC-AICtab$AIC[1])) # combined is the best AIC
  names(AIC) = c( "Model", "Covariates" , "df" , "AIC" , "deltaAIC")
  AIC


  # diagnostics
  mod=IGIYN_plus
  simulateResiduals(mod, plot= T)
  acf(residuals(mod))


  # statistics
  sum = summary(mod)
  stats =sum$coefficients$cond
  stats=apply (stats , 2 , function (x)  round (x,4))
  stats = cbind ( "foc.fami OES", stats)
  stats

  # save
  MOD3 = list ( AIC = AIC, stats = stats)
  save(MOD3, file= file.path(getwd() ,  "mod3.rda"))
} else {
  load ( file.path(getwd() ,  "mod3.rda"))
}

# stats
MOD3$AIC
MOD3$stats

# plot
boxplot( dat_den$foc.fami ~ dat_den$foc.oes)

###### MODEL 4 ######
runmodel4 = T
if(runmodel4){
  # models

  {
    IGITF_plus <- glmmTMB( riv.max.fami^0.2
                           ~
                             foc.oes+
                             igi_TF
                           + scale ( rain_md)
                           +  riv.max.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           ,data= dat_den
                           ,
                           beta_family(link = "logit")

    )
    print("1/8")

    IGITF_times <- glmmTMB( riv.max.fami^0.2
                            ~
                              foc.oes*
                              igi_TF
                            + scale ( rain_md)
                            +  riv.max.prop.coll.vs.hist
                            + comp
                            + (1| i.l.i_focal)
                            + ar1( datenF-1|acor)
                            , data= dat_den
                            ,
                            beta_family(link = "logit")
    )
    print("2/8")

    IGIFREQ_plus <- glmmTMB( riv.max.fami^0.2
                             ~
                               foc.oes+
                               igi_freq
                             + scale ( rain_md)
                             +  riv.max.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("3/8")

    IGIFREQ_times <- glmmTMB( riv.max.fami^0.2
                              ~
                                foc.oes*
                                igi_freq
                              + scale ( rain_md)
                              +  riv.max.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("4/8")

    nullmodel <- glmmTMB( riv.max.fami^0.2
                          ~
                            scale ( rain_md)
                          +  riv.max.prop.coll.vs.hist
                          + comp
                          + (1| i.l.i_focal)
                          + ar1( datenF-1|acor)
                          , data= dat_den
                          ,
                          beta_family(link = "logit")
    )
    print("5/8")

    IGIFREQ_only <- glmmTMB( riv.max.fami^0.2
                              ~
                                #foc.oes*
                                igi_freq
                              + scale ( rain_md)
                              +  riv.max.prop.coll.vs.hist
                              + comp
                              + (1| i.l.i_focal)
                              + ar1( datenF-1|acor)
                              , data= dat_den
                              ,
                              beta_family(link = "logit")
    )
    print("6/8")

    IGITF_only <- glmmTMB( riv.max.fami^0.2
                             ~
                               #foc.oes*
                               igi_TF
                             + scale ( rain_md)
                             +  riv.max.prop.coll.vs.hist
                             + comp
                             + (1| i.l.i_focal)
                             + ar1( datenF-1|acor)
                             , data= dat_den
                             ,
                             beta_family(link = "logit")
    )
    print("7/8")

   OES_only <- glmmTMB( riv.max.fami^0.2
                           ~
                             foc.oes

                           + scale ( rain_md)
                           +  riv.max.prop.coll.vs.hist
                           + comp
                           + (1| i.l.i_focal)
                           + ar1( datenF-1|acor)
                           , data= dat_den
                           ,
                           beta_family(link = "logit")
    )
    print("8/8")
  }

   # # bonus model
  IGITF_plusOES_plusBBY <- glmmTMB( riv.max.fami^0.2
                         ~
                           foc.oes+
                           igi_TF+
                           bby_TF
                         + scale ( rain_md)
                         +  riv.max.prop.coll.vs.hist
                         + comp
                         + (1| i.l.i_focal)
                         + ar1( datenF-1|acor)
                         ,data= dat_den
                         ,
                         beta_family(link = "logit")
  )
summary(IGITF_plusOES_plusBBY)
  # AIC
  AICtab = AIC(
    nullmodel,
               IGITF_times,
               IGIFREQ_times,
               IGITF_plus,
               IGIFREQ_plus,
               IGIFREQ_only,
               IGITF_only,
               OES_only)
  AICtab
  AICtab = AICtab[order(AICtab$AIC),]
  AIC= cbind("Rival familiarity OES", dimnames(AICtab)[[1]],AICtab,c(AICtab$AIC-AICtab$AIC[1])) # combined is the best AIC
  names(AIC) = c( "Model", "Covariates" , "df" , "AIC" , "deltaAIC")
  AIC

  #diagnostics
  mod = IGIFREQ_plus
  simulateResiduals(mod,plot=T)# not great
  hist(residuals  (mod)) # doesn't look too bad
  qqnorm(residuals(mod))# so try normal qq plot
  qqline(residuals(mod)) # looks fine

  # statistics
  sum = summary(mod)
  stats =sum$coefficients$cond
  stats=apply (stats , 2 , function (x)  round (x,4))
  stats = cbind ( "riv.fami OES", stats)


  # save
  MOD4 = list ( AIC = AIC, stats = stats)
  save(MOD4, file= file.path(getwd() ,  "mod4.rda"))
} else {
  load ( file.path(getwd() ,  "mod4.rda"))
}

# stats
MOD4$AIC
MOD4$stats

# plot
boxplot ( dat_den$riv.max.fami ~ dat_den$foc.oes)

### MODEL 5 ### Home range
runmodel5 = T
if( runmodel5){
hist ( dat_den$terr.size^0.4)
hist (  (log (  dat_den$terr.size)+5)^1.6 )

TERR<- glmmTMB( (log (  terr.size+1))^1.6
                ~
                  igi_TF
               # + bbyst_Y_N
               # + foc.oes
                + scale(rain_md)
                + comp
                + foc.prop.coll.vs.hist
                + (1| foc)
                + ar1( datenF-1|acor)
                , data= dat_den
                ,
                family= gaussian()
)


model_simres <- simulateResiduals(TERR, plot = T) # this is not great
hist( residuals(TERR)) #
ggplot ( dat_den , aes ( x = comp , y = terr.size, col = foc))+
  geom_point()+
  geom_smooth(method="lm") # broadly matches the model
sum = summary(TERR)
stats =sum$coefficients$cond
stats=apply (stats , 2 , function (x)  round (x,4))
stats = cbind ( "terr.size", stats)
stats
MOD5 = list ( stats = stats)
save(MOD5, file= file.path(getwd() ,  "mod5.rda"))
} else {
  load ( file.path(getwd() ,  "mod5.rda"))
}

#### MODEL 6 ### OVERLAP

runmodel6=T
if(runmodel6){
  # transformation
  #dat_den$trans_riv_max_overlap <-scale(dat_den$riv.max.overlap)
  dat_den$trans_riv_max_overlap <-dat_den$riv.max.overlap
  # dat_den$trans_riv_max_overlap <- dat_den$trans_riv_max_overlap+0+abs(min(na.omit(dat_den$trans_riv_max_overlap)))
  # dat_den$trans_riv_max_overlap <- dat_den$trans_riv_max_overlap / max( na.omit (dat_den$trans_riv_max_overlap))
  dat_den$trans_riv_max_overlap <- (dat_den$trans_riv_max_overlap+0.001)*0.999
  dat_den$trans_riv_max_overlap <- dat_den$trans_riv_max_overlap^0.7
  hist(dat_den$riv.max.overlap)
  hist(dat_den$trans_riv_max_overlap)
  range(na.omit(dat_den$trans_riv_max_overlap))

  # model
  OVERLAP <- glmmTMB( trans_riv_max_overlap
                      ~
                        riv.max.igi_freq+
                      bbyst_Y_N +
                        foc.oes
                      + scale(rain_md)
                      + foc.prop.coll.vs.hist
                      + riv.max.prop.coll.vs.hist
                      + (1| foc)
                      + (1|riv.max.id)
                      + ar1( datenF-1|acor)
                      , data= dat_den
                      ,
                      beta_family()
  )

  qqnorm(residuals(OVERLAP))
  qqline(residuals(OVERLAP))
  model_simres <- simulateResiduals(OVERLAP, plot = T) # this is not great, so errors are reported in main text
  hist( residuals(OVERLAP)) #
  sum = summary(OVERLAP)
  stats =sum$coefficients$cond
  stats=apply (stats , 2 , function (x)  round (x,4))
  stats = cbind ( "overlap", stats)
  stats
  MOD6 = list(stats = stats)
  save(MOD6, file= file.path(getwd() ,  "mod6.rda"))
} else {
  load ( file.path(getwd() ,  "mod6.rda"))
}

# plot
ggplot ( dat_den , aes( y = riv.max.overlap , x= riv.max.igi_freq, color = foc))+
  geom_point()+
  geom_smooth( method = "lm", se=F)





########## NOW COMBINE AND REPORT

bigstats=  rbind (
  MOD1$stats,
  MOD2$stats,
  MOD3$stats,
  MOD4$stats,
  MOD5$stats,
  MOD6$stats
)

bigAIC = rbind (
  MOD1$AIC,
  MOD2$AIC,
  MOD3$AIC,
  MOD4$AIC
)
write.csv( bigstats , file= file.path ( getwd(),  "statstable.csv"))
write.csv( bigAIC   , file= file.path ( getwd(), "AICtable.csv"))

