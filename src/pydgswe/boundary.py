class TransmissiveBoundary:
    def __call__(self, U_inside, U_inside_limit):
        return U_inside.const(), U_inside_limit
