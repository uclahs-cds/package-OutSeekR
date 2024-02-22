test_that(
    'least.significant.digit works', {
        expect_equal(
            least.significant.digit(
                x = c(0.2, 0.03, 0.004, 0.05)
                ),
            0.001
            );
        }
    );
