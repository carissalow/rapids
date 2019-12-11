Screen features
===============

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
