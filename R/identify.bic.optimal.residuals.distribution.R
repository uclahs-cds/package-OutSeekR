#' Identify optimal distribution of residuals
#'
#' Identify which of four distributions---normal, log-normal, exponential, or gamma---best fits the given vector of residuals according to BIC.
#'
#' @param x A numeric vector.
#' @return A numeric code representing which distribution optimally fits `x`.  Possible values are
#' * 1 = normal,
#' * 2 = log-normal,
#' * 3 = exponential, and
#' * 4 = gamma.
#' @export identify.bic.optimal.residuals.distribution
#' @examples
#' # Generate fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' identify.bic.optimal.residuals.distribution(
#'     x = x
#'     );
identify.bic.optimal.residuals.distribution <- function(x) {
    # Add a minimum value to ensure the values in `x` are strictly
    # positive.
    add.minimum.value <- least.significant.digit(x)
    if (min(x) < 0) {
        x.nozero <- x - min(x) + add.minimum.value;
        }
    else {
        x.nozero <- x + add.minimum.value;
        }

    # Fit a model of each candidate distribution to the data.
    glm.norm <- suppressWarnings(
        expr = gamlss::gamlss(
            formula = x.nozero ~ 1,
            family = gamlss.dist::NO(),
            control = gamlss::gamlss.control(
                trace = FALSE
                )
            )
        );
    glm.lnorm <- suppressWarnings(
        expr = gamlss::gamlss(
            formula = x.nozero ~ 1,
            family = gamlss.dist::LNO(),
            control = gamlss::gamlss.control(
                trace = FALSE
                )
            )
        );
    glm.exp <- suppressWarnings(
        expr = gamlss::gamlss(
            formula = x.nozero ~ 1,
            family = gamlss.dist::EXP(),
            control = gamlss::gamlss.control(
                trace = FALSE
                )
            )
        );
    glm.gamma <- suppressWarnings(
        expr = gamlss::gamlss(
            formula = x.nozero ~ 1,
            family = gamlss.dist::GA(),
            control = gamlss::gamlss.control(
                trace = FALSE
                )
            )
        );

    # Return the index of the distribution with the smallest BIC.
    which.min(c(glm.norm$sbc, glm.lnorm$sbc, glm.exp$sbc, glm.gamma$sbc));
    }
