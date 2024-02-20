import sqlite3

# Connect to the SQLite database
conn = sqlite3.connect('pathway-generation-sqlitedb')
cursor = conn.cursor()

# Insert a row into the table
cursor.execute('''
    INSERT INTO decomposed_uid_mapping (uid, component_id, input_or_output_uid, input_or_output_reactome_id, reactome_id)
    VALUES ('test_uid', 'test_component_id', 'test_input_output_uid', 1, 1)
''')
conn.commit()

# Print the contents of the table
print("Printing decomposed_uid_mapping table contents after insertion:")
cursor.execute('SELECT * FROM decomposed_uid_mapping')
for row in cursor.fetchall():
    print(row)

# Close the connection
conn.close()

