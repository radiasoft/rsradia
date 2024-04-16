
class Electromagnet():
    """A simple cylindrical ferromagnet wrapped in a solenoid"""
    
    def __init__(
            self, turns, hyst_model, M0=0., I0=0.
        ):
        """Args:
          * turns: Number of turns in the solenoid
          * hyst_model: Any rshyst hysteresis model
          * M0: Initial magnetization in the core (default 0.)
          * I0: Initial applied current (default 0.)
        """

        self.M = M0
        self.turns = turns
        self.model = hyst_model
        self.H = self.turns * I0

    def apply_current(self, I):
        """Applies a current to the solenoid, changing the magnet's field and magnetization"""
        
        H_path = [self.H, self.turns * I]
        self.H = H_path[1]
        _, M_path = self.model.path([H_path], self.M)
        self.M = M_path[-1]
    