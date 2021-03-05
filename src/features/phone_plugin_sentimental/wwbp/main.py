import pandas as pd
import numpy as np

def wwbp_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs): 

    sentiment_data = pd.read_csv(sensor_data_files["sensor_data"])
    sentiment_features = pd.DataFrame(columns=["local_segment"])

    # Get relevant information about provided features 
    app_included = "app" in provider["FEATURES"]

    if not sentiment_data.empty:

        sentiment_data = filter_data_by_segment(sentiment_data, time_segment)

        if not sentiment_data.empty:

            # Split the data into groups by time segment
            segments = sentiment_data.groupby("local_segment")
            time_segments = []

            # Get all tuples of app_category features to calculate
            if app_included:
                small_df = sentiment_data.drop(sentiment_data.columns.difference(['app_name', 'word_category']), axis=1)
                tuples = list(small_df.groupby(['app_name', 'word_category']).groups)
                categories = list(filter(lambda x : x[1] != 'total_words' and x[1]==x[1], tuples))

            # Get all word categories features to calculate 
            else:
                categories = sentiment_data.word_category.unique().tolist()
                categories.remove('total_words')

            # Aggregate each segment data into a single instance 
            for _, segment_df in segments:
                instance = process_local_segment(segment_df, categories, app_included)
                time_segments.append(instance)

            # Combine the data into a final dataframe
            sentiment_features = pd.concat(time_segments)

    return sentiment_features


# Handles local segment logic and processing
def process_local_segment(df, categories, app_included):

    if not app_included:
        # Simply add the relevant features
        features_df = insert_features(df, categories, False)

    else:
        # Get the features for each app-category combination 
        app_groups = df.groupby(['app_name'])
        processed = []
        for _, group in app_groups:
            p = insert_features(group, categories, True)
            processed.append(p)

        # Combine the data into one instance 
        features_df = pd.concat(processed).groupby(['local_segment'], as_index = False).sum()
    
    # Add the device_id column and return the data
    features_df['device_id'] = df['device_id'].values[0]

    return features_df


# Calculates features in a particular local segment 
def insert_features(df, categories, app_included=False):

    app = df['app_name'].values[0]

    # Map each word_category to its score 
    category_to_score = {}
    totals_df = pd.DataFrame()
    totals_df['score'] = df.groupby(["word_category"])['double_sentiment_score'].sum()
    for index, row in totals_df.iterrows():                
        category_to_score[index] = row['score']

    # Get the total number of words in the time segment
    total_words = category_to_score['total_words']
    category_to_score.pop('total_words')

    # Populate data with the available scores otherwise fill in 0
    data = {}
    for c in categories:

        # c is a tuple (app, word_category)
        if app_included:
            tuple_app = c[0]
            tuple_cat = c[1]

            # Calculate the normalized score if c present 
            feature = tuple_cat + "_" + tuple_app
            if tuple_app == app and tuple_cat in category_to_score:
                data[feature] = category_to_score[tuple_cat] / total_words
            else:
                data[feature] = 0

        # c is just a word_category
        else: 
            # Calculate the normalized score if c is present 
            feature = c 
            if c in category_to_score:
                data[feature] = category_to_score[c] / total_words
            else:
                data[feature] = 0

    # Create a dataframe from the data
    data['local_segment'] = df['local_segment'].values[0]
    processed_df = pd.DataFrame([data])

    return processed_df

