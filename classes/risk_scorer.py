from classes.logger import Logger


class RiskScorer:
    def __init__(self, weights=None):
        self.logger = Logger()
        self.weights = weights or {}
        self.logger.info("Risk scorer initialized.")

    def set_weights(self, weights):
        if not isinstance(weights, dict):
            self.logger.failure("Risk weights must be provided as a dictionary.")
            raise TypeError("Weights must be a dictionary.")

        if not weights:
            self.logger.failure("Risk weights cannot be empty.")
            raise ValueError("Weights cannot be empty.")

        for feature, weight in weights.items():
            if not isinstance(weight, (int, float)):
                self.logger.failure(
                    f"Risk weight for {feature} is not numeric."
                )
                raise TypeError("Risk weights must be numeric.")

            if weight < 0:
                self.logger.failure(
                    f"Risk weight for {feature} is negative."
                )
                raise ValueError("Risk weights cannot be negative.")

        if sum(weights.values()) == 0:
            self.logger.failure("Risk weights cannot all be zero.")
            raise ValueError("Risk weights cannot all be zero.")

        self.weights = weights.copy()

        self.logger.info(
            f"Risk weights configured for {len(self.weights)} features."
        )

    def get_weights(self):
        return self.weights.copy()

    def normalize_weights(self):
        if not self.weights:
            self.logger.failure("Cannot normalize empty risk weights.")
            raise ValueError("Risk weights cannot be empty.")

        total = sum(self.weights.values())

        if total == 0:
            self.logger.failure("Cannot normalize zero-sum risk weights.")
            raise ValueError("Risk weights cannot all be zero.")

        normalized = {
            feature: weight / total
            for feature, weight in self.weights.items()
        }

        self.weights = normalized

        self.logger.info("Risk weights normalized.")

        return normalized

    def score_features(self, features):
        if not isinstance(features, dict):
            self.logger.failure("Risk features must be provided as a dictionary.")
            raise TypeError("Features must be a dictionary.")

        if not self.weights:
            self.logger.failure("Risk scoring attempted without configured weights.")
            raise ValueError("Risk weights must be configured before scoring.")

        missing_features = [
            feature for feature in self.weights
            if feature not in features
        ]

        if missing_features:
            self.logger.failure(
                f"Risk scoring failed: missing features {missing_features}."
            )
            raise KeyError(
                f"Missing features required for scoring: {missing_features}."
            )

        normalized_weights = self.weights.copy()
        total_weight = sum(normalized_weights.values())

        if total_weight == 0:
            self.logger.failure("Risk scoring failed: total weight is zero.")
            raise ValueError("Risk weights cannot all be zero.")

        if total_weight != 1:
            normalized_weights = {
                feature: weight / total_weight
                for feature, weight in normalized_weights.items()
            }

        contributions = {}

        for feature, weight in normalized_weights.items():
            value = features[feature]

            if not isinstance(value, (int, float)):
                self.logger.failure(
                    f"Risk feature {feature} is not numeric."
                )
                raise TypeError("Risk features must be numeric.")

            value = min(max(value, 0.0), 1.0)

            contributions[feature] = value * weight

        score = sum(contributions.values())
        score = min(max(score, 0.0), 1.0)

        self.logger.info(
            f"Risk score calculated: {score:.4f}."
        )

        return {
            "score": score,
            "contributions": contributions
        }

    def score_regions(self, regional_features):
        if not isinstance(regional_features, list):
            self.logger.failure(
                "Regional risk scoring requires a list of regions."
            )
            raise TypeError("Regional features must be a list.")

        scored_regions = []

        for region in regional_features:
            if not isinstance(region, dict):
                self.logger.failure(
                    "Regional risk scoring encountered a non-dictionary region."
                )
                raise TypeError("Each region must be a dictionary.")

            features = region.get("features")

            if not isinstance(features, dict):
                self.logger.failure(
                    "Regional risk scoring encountered a region without features."
                )
                raise KeyError("Each region must contain a features dictionary.")

            result = self.score_features(features)

            scored_region = region.copy()
            scored_region["risk_score"] = result["score"]
            scored_region["contributions"] = result["contributions"]

            scored_regions.append(scored_region)

        self.logger.info(
            f"Risk scores calculated for {len(scored_regions)} regions."
        )

        return scored_regions