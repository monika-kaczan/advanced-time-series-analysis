######### Biblioteki #########

library(dplyr)
library(xts)
library(fBasics)
library(tseries)
library(car)
library(FinTS)
library(fGarch)
library(rugarch)



######### Wczytanie danych #########

# Roczne dane
DASH <- read.csv("DASH.csv",
                   header = TRUE,
                   sep = ";",
                   dec = ".",
                   stringsAsFactors = F)
KCS <- read.csv("KCS.csv",
                 header = TRUE,
                 sep = ";",
                 dec = ".",
                 stringsAsFactors = F)
MKR <- read.csv("MKR.csv",
                header = TRUE,
                sep = ";",
                dec = ".",
                stringsAsFactors = F)
MATIC <- read.csv("MATIC.csv",
                  header = TRUE,
                  sep = ";",
                  dec = ".",
                  stringsAsFactors = F)


marketCap <- data.frame(DASH$marketCap, KCS$marketCap, MATIC$marketCap, MKR$marketCap)

portofolio <- ( DASH$marketCap / rowSums(marketCap) * DASH$close ) +
  ( KCS$marketCap / rowSums(marketCap) * KCS$close ) +
  ( MATIC$marketCap / rowSums(marketCap) * MATIC$close ) +
  ( MKR$marketCap / rowSums(marketCap) * MKR$close )

plot(as.Date(DASH$timestamp), MATIC$marketCap/1000000, type = "l", col = "red", xlab = "Date", ylab = "Market capitalization [mln $]")
lines(as.Date(DASH$timestamp), KCS$marketCap/1000000, type = "l", col = "blue")
lines(as.Date(DASH$timestamp), DASH$marketCap/1000000, type = "l", col = "purple")
lines(as.Date(DASH$timestamp), MKR$marketCap/1000000, type = "l", col = "green")
mtext("% Share", side = 4, line = 3)
legend("topright", legend=c("MATIC", "KCS", "DASH", "MKR"),
       col=c("red", "blue", "purple", "green"), lty = 1, cex = 0.8)

crypto <- data.frame(portofolio, as.Date(DASH$timestamp))
colnames(crypto) <- c("crypto", "date")

# Jako okres przyjmujemy ostatnie 2,5 roku, dokładniej od stycznia 2021. Pozwala to na uzyskanie
# odpowiedniej liczby obserwacji, a jednocześnie zachowanie analizy w miarę aktualnej.




######### Analiza wstepna #########



# Obliczamy ciagle zwroty
crypto$r <- diff.xts(log(crypto$crypto))


# Wykres zwrotow i notowan

par(mfrow = c(2, 1))
plot(crypto$date, crypto$r,
     type = "l", col="red", lwd = 1, xlab = "Data", ylab = "",
     main = "Zwroty portfela")
plot(crypto$date, crypto$crypto,
     type = "l", col = "black", lwd = 1, xlab = "Data", ylab = "$",
     main = "Notowania portfela")
par(mfrow = c(1, 1))



# Czy zwroty pochodzą z rozkładu normalnego?

hist(crypto$r, prob = T, breaks = 100, main = "Histogram zwrotów")
curve(dnorm(x,
            mean = mean(crypto$r, na.rm = T),
            sd = sd(crypto$r, na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)

# Widzimy, że mamy relatywnie dużo wysokich wartości, jednak rozkład nie ma
# np. wysokiego szczytu funkcji gęstości charakterytycznego dla rozkładu leptokurtycznego.

# Porównajmy rozkład zwrotów i normalny za pomocą QQ plot.

qqnorm(crypto$r)
qqline(crypto$r, col = 2)

# Widzimy duże odchodzenie od roz. normalnego w ogonach.

# Możemy też formalnie sprawdzić to np. testem Jacque-Bera.

jarque.bera.test(na.omit(crypto$r))

# H0 o normalności silnie odrzucana.


# Sprawdźmy również autokorelację za pomocą testu DW

set.seed(410998) # Reproduktowalność wyników testu DW
durbinWatsonTest(lm(formula = crypto$r ~ 1),
                 max.lag = 10)

# Na poziomie istoności 5% mamy autokorelację dla 2 i 8 opóźnienia



# Wykres ACF zwrotow

acf(crypto$r, lag.max = 100, na.action = na.pass,
    col = "darkblue", lwd = 7,
    main = "Wykres ACF zwrotow crypto")

# Przybliżmy
acf(crypto$r, lag.max = 100, na.action = na.pass,
    ylim = c(-0.2, 0.2),
    col = "darkblue", lwd = 7,
    main = "Wykres ACF zwrotow crypto")

# Widzimy kilka istotnych wartośći, co może wskazywać na pewne problemy z autokorelacją.


# Wykres ACF kwadratow zwrotow

acf((crypto$r^2), lag.max = 100, na.action = na.pass,
    col = "darkblue", lwd = 7,
    main = "Wykres ACF kwadratow zwrotow crypto")

# Widzimy kilka istotnych wartości, jednak niedużo, brak zależności korelacyjnych wyżnych
# rzędów i nie wygasają one. -> brak wskazań że mamy efekty ARCH
# [? nie jestem tego pewna właśnie to nie wychodzi na wykresie :(]


# Ponieważ na wykresie nie wychodzi, sprawdźmy występowanie efektów ARCH wśród zwrotów
# formalnym testem.

ArchTest(crypto$r, lags = 5)
ArchTest(crypto$r, lags = 10)

# W obu przypadkach odrzucamy H0 o braku efektów ARCH -> możemy modelować za pomocą modeli GARCH.





######### Modelowanie #########

# 1.
# ARCH(1)
k.arch1 <- garchFit(formula = ~ garch(1, 0),
                    data = na.omit(crypto$r),
                    cond.dist = "norm", # rozkład warunkowy reszt
                    trace = FALSE) # jeśli nie chcemy oglądać szczegółów poszczególnych iteracji

# podsumowanie wyników i kilka testów diagnostycznych
summary(k.arch1)

# stała w równaniu średniej jest nieistotna, więc możemy ją wykluczyć z modelu
k.arch1 <- garchFit(formula = ~ garch(1, 0),
                    data = na.omit(crypto$r),
                    include.mean = F,
                    cond.dist = "norm",
                    trace = FALSE)

# podsumowanie wyników i kilka testów diagnostycznych
summary(k.arch1)

# Wartości kryteriów informacyjnych są nieznacznie niższe, więc poprawiliśmy model.

# Na podstawie testu JB Widzimy, że reszty nie mają rozkładu normalnego.

# Na podstawie testu LM ARCh widzimy, że nie udało nam się wyeliminować efektów ARCH
# (odrzucamy H0 o braku efektów ARCH)


# 2.
# ARCH(3)
k.arch3 <- garchFit(formula = ~ garch(3, 0),
                    data = na.omit(crypto$r),
                    include.mean = F,
                    cond.dist = "norm",
                    trace = FALSE)
summary(k.arch3)

# Wszystkie parametry są statystycznie istotne. Udało nam się też usunąć efekty ARCH i poprawiliśmy
# wartości kryteriów informacyjnych.

# Sprójrzmy na wykresy ACF dla standaryzowanych reszt
# i ich kwadratów (10 i 11)
# plot(k.arch3)
# Wszystkie statystycznie nieistotne więc super.


# 4.
# GARCH(1,1)
k.garch11 <- garchFit(formula = ~ garch(1, 1),
                      data = na.omit(crypto$r),
                      include.mean = F,
                      cond.dist = "norm",
                      trace = FALSE)
summary(k.garch11)

# Poprawiliśmy wartości kryteriów informacyjnych, ale nieznacznie. Ponownie mamy
# wyeliminowane efekty ARCH.

# Sprójrzmy na wykresy ACF dla standaryzowanych reszt
# i ich kwadratów (10 i 11)
plot(k.garch11)

# Tutaj zmienialiśmy wartości parametrów i sprawdzaliśmy wyniki
k.garch12 <- garchFit(formula = ~ garch(3, 1),
                      data = na.omit(crypto$r),
                      cond.dist = "norm",
                      trace = FALSE)
summary(k.garch12)




######### Oszacowania funkcji warunkowej wariancji #########

k.garch11@fit$par

var_uncond <- k.garch11@fit$par[1] / (1 - k.garch11@fit$par[2]
                                      - k.garch11@fit$par[3])
names(var_uncond) <- "unconditional variance"
var_uncond

# szacujemy prognoze warunkowej wariancji na nastepne 100 okresow
# wyniki modelu w k.garch11
k.fore100 <- predict(k.garch11, n.ahead = 100)
head(k.fore100)

# ostatnia kolumna (trzecia) zawiera prognozy warunkowego odchylenia standardowego
# musimy podniesc je do kwadratu aby otrzymac prognozy warunkowej wariancji

# wykres oszacowan i prognoz warunkowej wariancji w dlugim okresie (bierzemy odchylenie stand to kwadratu czyliw wariancę)
plot(k.fore100[, 3] ^ 2, type = "l")
# dodajmy jeszcze poziom wariancji bezwarunkowej
abline(h = var_uncond, col = "red", lty = 2)
title(main = "Warunkowa i bezwarunkowa wariancja zwrotów WIG20")
# (Widzimy tylko dolną linię z wykresu na s. 17 w slajdach) [# nie wiem co mamy jeszcze widzieć?]
# plot(k.fore100[, 3] ^ 2, type = "l", ylim = c(0.004, 0.02))
# abline(h = var_uncond, col = "red", lty = 2)

# powtorzmy procedure dla okresu skladajacego sie z 500 interwalow
k.fore500 <- predict(k.garch11, n.ahead = 500)
plot(k.fore500[, 3] ^ 2, type = "l")
abline(h = var_uncond, col = "red", lty = 2)
title(main = "Warunkowa i bezwarunkowa wariancja zwrotów WIG20")

# Mozemy zauwazyc, ze prognozy warunkowej wariancji zbiegaja w dlugim okresie
# do poziomu wariancji bezwarunkowej.


## Garch-T(1,1)

spec3 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                    distribution.model = "std")
portfel2.garcht11 <- ugarchfit(spec=spec3, data=na.omit(crypto$r))
portfel2.garcht11


plot(portfel2.garcht11, which = 10)


# t-GARCH(1,1)

spec2 <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(2, 1), submodel = "TGARCH"),
                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                    distribution.model = "norm")
portfel2.tgarch11 <- ugarchfit(spec = spec2, data=na.omit(crypto$r))
portfel2.tgarch11


plot(portfel2.tgarch11, which = 10)

## gjrGARCH


spec2 <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1), submodel = "TGARCH"),
                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                    distribution.model = "norm")
portfel2.tgarch11 <- ugarchfit(spec = spec2, data=na.omit(crypto$r))
portfel2.tgarch11


plot(portfel2.tgarch11, which = 10)


## Wartość narażona na ryzyko

data_in_sample <- crypto[1:700,]

# standaryzacja zwrotów
data_in_sample$rstd <- (data_in_sample$r - mean(data_in_sample$r, na.rm=T)) /
  sd(data_in_sample$r ,na.rm = T)

q01 <- quantile(data_in_sample$rstd, 0.01, na.rm = T)


## Wyliczanie modeli

## ARCH(1)

spec0 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = T),
                    distribution.model = "norm")
portfel2.arch1 <- ugarchfit(spec = spec0, data = na.omit(data_in_sample$r))
portfel2.arch1

## GARCH(1,1)

spec0 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = T),
                    distribution.model = "norm")
portfel2.garch1 <- ugarchfit(spec = spec0, data = na.omit(data_in_sample$r))
portfel2.garch1

## GARCH-t(1,1)

spec3 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                    distribution.model = "std")
portfel2.garcht11 <- ugarchfit(spec=spec3, data=na.omit(data_in_sample$r))
portfel2.garcht11

# t-GARCH(1,1)

spec2 <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(2, 1), submodel = "TGARCH"),
                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                    distribution.model = "norm")
portfel2.tgarch11 <- ugarchfit(spec = spec2, data=na.omit(data_in_sample$r))
portfel2.tgarch11


## gjrGARCH


spec2 <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1), submodel = "TGARCH"),
                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                    distribution.model = "norm")
portfel2.gjrgarch11 <- ugarchfit(spec = spec2, data=na.omit(data_in_sample$r))
portfel2.gjrgarch11


## Liczenie
data_in_sample <- data_in_sample[2:700,]

data_in_sample$VaR <- q01 * portfel2.arch1@fit$sigma
data_in_sample$VaR2 <- q01 * portfel2.garch1@fit$sigma
data_in_sample$VaR3 <- q01 * portfel2.garcht11@fit$sigma
data_in_sample$VaR4 <- q01 * portfel2.tgarch11@fit$sigma
data_in_sample$VaR5 <- q01 * portfel2.gjrgarch11@fit$sigma


# wyres zwrotóW vs. VaR
plot(data_in_sample$date, data_in_sample$r, col = "red", lwd = 1, type = 'l',
     ylim = c(-0.5, 0.5), main = "ARCH(1)")
lines(data_in_sample$date, data_in_sample$VaR, type = 'l', col = "green")

plot(data_in_sample$date, data_in_sample$r, col = "red", lwd = 1, type = 'l',
     ylim = c(-0.5, 0.5), main = "GARCH(1,1)")
lines(data_in_sample$date, data_in_sample$VaR2, type = 'l', col = "green")

plot(data_in_sample$date, data_in_sample$r, col = "red", lwd = 1, type = 'l',
     ylim = c(-0.5, 0.5), main = "GARCH-t(1,1)")
lines(data_in_sample$date, data_in_sample$VaR3, type = 'l', col = "green")

plot(data_in_sample$date, data_in_sample$r, col = "red", lwd = 1, type = 'l',
     ylim = c(-0.5, 0.5), main = "tGARCH(2,1)")
lines(data_in_sample$date, data_in_sample$VaR4, type = 'l', col = "green")

plot(data_in_sample$date, data_in_sample$r, col = "red", lwd = 1, type = 'l',
     ylim = c(-0.5, 0.5), main = "gjr-GARCH(1,1)")
lines(data_in_sample$date, data_in_sample$VaR5, type = 'l', col = "green")


# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
sum(data_in_sample$r < data_in_sample$VaR) / length(data_in_sample$VaR)

sum(data_in_sample$r < data_in_sample$VaR2) / length(data_in_sample$VaR)

sum(data_in_sample$r < data_in_sample$VaR3) / length(data_in_sample$VaR)

sum(data_in_sample$r < data_in_sample$VaR4) / length(data_in_sample$VaR)

sum(data_in_sample$r < data_in_sample$VaR5) / length(data_in_sample$VaR)


## Out-of-sample

start  <- 632
finish <- 731

Portfel4 <- crypto[632:731, ]

VaR <- rep(NA, times = finish - start + 1)
VaR2 <- rep(NA, times = finish - start + 1)
VaR3 <- rep(NA, times = finish - start +1)

crypto$obs <- 1:length(crypto$r)

spec <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(2, 1), submodel = "TGARCH"),
                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                   distribution.model = "norm")


for (k in start:finish) {
  tmp.data <- crypto[crypto$obs <= (k - 1), ]
  tmp.data <- tmp.data[as.Date("2022-11-16") <= tmp.data$date, ]
  tmp.data$rstd <- (tmp.data$r - mean(tmp.data$r, na.rm = T)) /
    sd(tmp.data$r, na.rm = T)
  q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
  tmp.egarch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$r))
  sigma.forecast  <- ugarchforecast(tmp.egarch11, n.ahead = 1)
  sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
  VaR[k - start + 1] <- q01 * sigma.forecast2
}

Portfel4$VaR <- VaR


plot(Portfel4$date, Portfel4$r, col = "red", lwd = 1, type = 'l',
     ylim = c(-0.25, 0.25), main ="t-GARCH(2,1)")
abline(h = 0, lty = 2)
lines(Portfel4$date, Portfel4$VaR, type = 'l', col = "green")

# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
sum(Portfel4$r < Portfel4$VaR) / length(Portfel4$VaR)