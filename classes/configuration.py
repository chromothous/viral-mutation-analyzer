from classes.logger import Logger


class Configuration:
    def __init__(self):
        self.logger = Logger()
        self.settings = {
            "window_size": 100,
            "kmer_size": 3,
            "repeat_length": 3,
            "motifs": {}
        }
        self.logger.info("Configuration initialized.")

    def set(self, name, value):
        if name not in self.settings:
            self.logger.failure(
                f"Configuration setting does not exist: {name}."
            )
            raise KeyError(f"Unknown configuration setting: {name}.")

        self.settings[name] = value

        self.logger.info(
            f"Configuration setting updated: {name}."
        )

    def get(self, name):
        if name not in self.settings:
            self.logger.failure(
                f"Configuration setting does not exist: {name}."
            )
            raise KeyError(f"Unknown configuration setting: {name}.")

        return self.settings[name]

    def get_all(self):
        return self.settings.copy()

    def reset(self):
        self.settings = {
            "window_size": 100,
            "kmer_size": 3,
            "repeat_length": 3,
            "motifs": {}
        }
        self.logger.info("Configuration reset to defaults.")