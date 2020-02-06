Extracted Features
==================

Accelerometer
-------------

Available epochs: daily, morning, afternoon, evening, and night

- Max magnitude: maximum magnitude of acceleration (:math:`\|acceleration\| = \sqrt{x^2 + y^2 + z^2}`)
- Min magnitude: minimum magnitude of acceleration
- Avg magnitude: average magnitude of acceleration
- Median magnitude: median magnitude of acceleration
- Std magnitude: standard deviation of acceleration
- Ratio exertional activity episodes: ratio of exertional activity time periods to total time periods
- Sum exertional activity episodes: total minutes of performing exertional activity during the epoch
- Longest exertional activity episode: longest episode of performing exertional activity
- Longest non-exertional activity episode: longest episode of performing non-exertional activity
- Count exertional activity episodes: count of the episods of performing exertional activity
- Count non-exertional activity episodes: count of the episodes of performing non-exertional activity

Battery
--------

Available epochs: daily, morning, afternoon, evening, and night

-	Count discharge: number of battery discharging episodes
-	Sum duration discharge: total duration of all discharging episodes (time the phone was discharging)
-	Average consumption rate: average of the ratios between discharging episodes’ battery delta and duration
-	Max consumption rate: max of the ratios between discharging episodes’ battery delta and duration
-	Count charge: number of battery charging episodes
-	Sum duration charge: total duration of all charging episodes (time the phone was charging)

Bluetooth
---------

Available epochs: daily, morning, afternoon, evening, and night

-	Count of scans (a scan is a row containing a single Bluetooth device detected by Aware)
-	Unique devices (number of unique devices identified by their hardware address -bt_address field)
-	Count of scans of the most unique device across each participant’s dataset 

Calls
-----

Available epochs: daily, morning, afternoon, evening, and night

-	Outgoing: count, count of distinct contacts, mean duration, sum duration, min duration, max duration, std duration, mode duration, entropy duration, time of first call (hours), time of last call (hours), count of most frequent contact.
-	Received: count, count of distinct contacts, mean duration, sum duration, min duration, max duration, std duration, mode duration, entropy duration, time of first call (hours), time of last call (hours), count of most frequent contact.
-	Missed: count, distinct contacts, time of first call (hours), time of last call (hours), count of most frequent contact.

Google Activity Recognition
---------------------------

Available epochs: daily, morning, afternoon, evening, and night

-	Count (number of rows)
-	Most common activity
-	Number of unique activities
-	Activity change count (any transition between two different activities, sitting to running for example)
-	Sum stationary: total duration of episodes of still and tilting (phone) activities
-	Sum mobile: total duration of episodes of on foot, running, and on bicycle activities
-	Sum vehicle: total duration of episodes of on vehicle activity

Light
-----

Available epochs: daily, morning, afternoon, evening, and night

- Count (number of rows)
- Max lux: maximum ambient luminance in lux units
- Min lux: minimum ambient luminance in lux units
- Avg lux: average ambient luminance in lux units
- median lux: median ambient luminance in lux units
- Std lux: standard deviation of ambient luminance in lux units

Location (Barnett’s) Fetures
-----------------------------

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

Screen
------

Available epochs: daily, morning, afternoon, evening, and night

Notes. An unlock episode is considered as the time between an unlock event and a lock event. iOS recorded these episodes reliable (albeit duplicated lock events within milliseconds from each other). However, in Android there are multiple consecutive unlock/lock events so we keep the closest pair. This happens because ACTION_SCREEN_OFF and ON are "sent when the device becomes non-interactive which may have nothing to do with the screen turning off" see this link

-	Count on: count of screen on events (only available for Android)
-	Count unlock: count of screen unlock events
-	Diff count on off: For debug purposes, on and off events should come in pairs, difference should be close to zero then.
-	Diff count unlock lock, For debug purposes, unlock and lock events should come in pairs, difference should be close to zero then.
-	Sum duration unlock: sum duration of unlock episodes 
-	Max duration unlock: maximum duration of unlock episodes
-	Min duration unlock: minimum duration of unlock episodes
-	Average duration unlock: average duration of unlock episodes
-	Std duration unlock: standard deviation of the duration of unlock episodes

SMS
---

Available epochs: daily, morning, afternoon, evening, and night

-	Sent: count, distinct contacts, time first sms, time last sms, count most frequent contact
-	Received: count, distinct contacts, time first sms, time last sms, count most frequent contact
