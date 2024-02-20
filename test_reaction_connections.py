import os
import sqlite3

print("T E S T")

try:
    # Connect to the database
    conn = sqlite3.connect('pathway_data.db')
    c = conn.cursor()

    # Check if the reaction_connections table exists
    c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='reaction_connections'")
    result = c.fetchone()

    if result:
        print("The reaction_connections table exists.")
        # Execute a query to select all rows from the reaction_connections table
        c.execute('SELECT * FROM reaction_connections')

        # Fetch all rows
        rows = c.fetchall()

        # Print the rows
        for row in rows:
            print(row)

        # Execute a query to count the number of rows in the reaction_connections table
        c.execute('SELECT COUNT(*) FROM reaction_connections')

        # Fetch the count
        count = c.fetchone()[0]

        # Print the count
        print("Number of rows in reaction_connections table:", count)

    else:
        print("The reaction_connections table does not exist.")
       
    # Close the connection
    conn.close()
except sqlite3.Error as e:
    print("SQLite error:", e)
except Exception as e:
    print("An error occurred:", e)

# Print the current working directory to check the script's location
print("Current directory:", os.getcwd())


