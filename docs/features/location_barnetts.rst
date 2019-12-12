Location (Barnett’s) Fetures
=============================

Available epochs: daily

Barnett’s location features are based on the concept of flights and pauses. GPS coordinates are converted into a sequence of flights (straight line movements) and pauses (time spent stationary). Data is imputed before metrics are computed (https://arxiv.org/abs/1606.06328)

-	Time at home. Time spent at home in minutes. Home is the most visited significant location between 8 pm and 8 am including any pauses within a 200-meter radius.
-	Max home distance. Maximum distance from home in meters.
-	Pause probability. The fraction of a day spent in a pause (as opposed to a flight)
-	Circadian routine. A continuous metric that can take any value between 0 and 1, where 0 represents a daily routine completely different from any other sensed days and 1 a routine the same as every other sensed day.
-	Wkn circadian routine. Same as Circadian routine but computed separately for weekends and weekdays.
-	Distance travelled. Total distance travelled over a day.
-	Radius of Gyration (RoG). It is a measure in meters of the area covered by a person over a day. A centroid is calculated for all the places (pauses) visited during a day and a weighted distance between all places and the centroid is computed. The weights are proportional to the time spent in each place.
-	Maximum diameter. Largest distance in meters between any two pauses.
-	Avg flight duration. Mean duration of all flights.
-	Avg flight length. Mean length of all flights
-	Std flight duration. The standard deviation of the duration of all flights.
-	Std flight length. The standard deviation of the length of all flights.
-	Significant locations. The number of significant locations visited during the day. Significant locations are computed using k-means clustering over pauses found in the whole monitoring period. The number of clusters is found iterating from 1 to 200 stopping until the centroids of two significant locations are within 400 meters of one another.
-	Significant location entropy. Entropy measurement based on the proportion of time spent at each significant location visited during a day.
