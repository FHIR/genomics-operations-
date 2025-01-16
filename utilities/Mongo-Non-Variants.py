import os
import json
from pymongo import MongoClient
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# MongoDB connection details
MONGO_URI = os.getenv("MONGO_URI", "mongodb://localhost:27017")
DATABASE_NAME = "your_database_name"  # Replace with your database name

# Folder containing JSON files
OUTPUT_FOLDER = "NonVariantFolder"  # Replace with the actual path

# Ensure the output folder exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)


def push_to_mongo(database):
    """Push JSON data from the folder to MongoDB collections."""
    bed_collection = database.Bed
    non_variant_collection = database.NonVariant

    # Loop through all JSON files in the folder
    for file_name in os.listdir(OUTPUT_FOLDER):
        if file_name.endswith(".json"):
            file_path = os.path.join(OUTPUT_FOLDER, file_name)
            try:
                # Read JSON data from file
                with open(file_path, "r") as file:
                    data = json.load(file)

                # Insert into MongoDB
                if isinstance(data, list):  # Ensure data is list-like for bulk insertion
                    if "bed" in file_name.lower():
                        bed_collection.insert_many(data)
                        print(f"Inserted {len(data)} records into Bed collection from {file_name}.")
                    elif "non_variant" in file_name.lower():
                        non_variant_collection.insert_many(data)
                        print(f"Inserted {len(data)} records into NonVariant collection from {file_name}.")
                else:
                    print(f"Skipped {file_name}: Data is not a list and cannot be inserted.")
            except Exception as e:
                print(f"Error processing {file_name} for MongoDB: {e}")


def main():
    # Connect to MongoDB
    client = MongoClient(MONGO_URI)
    database = client[DATABASE_NAME]

    try:
        # Push data to MongoDB
        push_to_mongo(database)
    except Exception as e:
        print(f"An error occurred while pushing to MongoDB: {e}")
    finally:
        # Close MongoDB connection
        client.close()
        print("MongoDB connection closed.")


if __name__ == "__main__":
    main()
