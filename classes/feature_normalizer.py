from classes.logger import Logger


class FeatureNormalizer:
    def __init__(self):
        self.logger = Logger()
        self.logger.info("Feature normalizer initialized.")

    def normalize_value(self, value, minimum, maximum, invert=False):
        if not isinstance(value, (int, float)):
            self.logger.failure("Normalization failed: value is not numeric.")
            raise TypeError("Value must be numeric.")

        if not isinstance(minimum, (int, float)):
            self.logger.failure("Normalization failed: minimum is not numeric.")
            raise TypeError("Minimum must be numeric.")

        if not isinstance(maximum, (int, float)):
            self.logger.failure("Normalization failed: maximum is not numeric.")
            raise TypeError("Maximum must be numeric.")

        if minimum > maximum:
            self.logger.failure("Normalization failed: minimum exceeds maximum.")
            raise ValueError("Minimum cannot exceed maximum.")

        if minimum == maximum:
            self.logger.warning(
                "Normalization range contains identical minimum and maximum values."
            )
            return 0.0

        normalized = (value - minimum) / (maximum - minimum)
        normalized = min(max(normalized, 0.0), 1.0)

        if invert:
            normalized = 1.0 - normalized

        return normalized

    def normalize_features(self, features, ranges, inverted_features=None):
        if not isinstance(features, dict):
            self.logger.failure("Feature normalization failed: features are not a dictionary.")
            raise TypeError("Features must be provided as a dictionary.")

        if not isinstance(ranges, dict):
            self.logger.failure("Feature normalization failed: ranges are not a dictionary.")
            raise TypeError("Ranges must be provided as a dictionary.")

        if inverted_features is None:
            inverted_features = set()

        normalized = {}

        for feature, value in features.items():
            if feature not in ranges:
                self.logger.failure(
                    f"Feature normalization failed: range missing for {feature}."
                )
                raise KeyError(f"Normalization range missing for feature: {feature}.")

            minimum = ranges[feature]["minimum"]
            maximum = ranges[feature]["maximum"]
            invert = feature in inverted_features

            normalized[feature] = self.normalize_value(
                value,
                minimum,
                maximum,
                invert
            )

        self.logger.info(
            f"Normalized {len(normalized)} features."
        )

        return normalized

    def normalize_regional_features(
        self,
        regional_features,
        ranges,
        inverted_features=None
    ):
        if not isinstance(regional_features, list):
            self.logger.failure(
                "Regional feature normalization failed: input is not a list."
            )
            raise TypeError("Regional features must be provided as a list.")

        normalized_regions = []

        for region in regional_features:
            if not isinstance(region, dict):
                self.logger.failure(
                    "Regional feature normalization failed: region is not a dictionary."
                )
                raise TypeError("Each region must be a dictionary.")

            normalized_region = region.copy()
            normalized_region["features"] = self.normalize_features(
                region.get("features", {}),
                ranges,
                inverted_features
            )
            normalized_regions.append(normalized_region)

        self.logger.info(
            f"Normalized features for {len(normalized_regions)} regions."
        )

        return normalized_regions