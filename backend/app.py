from flask import Flask, request, jsonify, send_file, send_from_directory
import os
from werkzeug.utils import secure_filename
app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
RESULT_FOLDER = 'results'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULT_FOLDER'] = RESULT_FOLDER

@app.route('/')
def home():
    return send_from_directory('../frontend', 'index.html')

@app.route('/<path:filename>')
def static_files(filename):
    return send_from_directory('../frontend', filename)

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400
    file = request.files['file']
    filename = secure_filename(file.filename)
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(file_path)
    # Run the model here (placeholder)
    result_path = temp_model(file_path)
    return jsonify({'result_file': result_path})


#TODO: Implement the data annotation tool and predictive model here.
def temp_model(input_file):
    # Just copy the input file to the results folder
    os.makedirs(RESULT_FOLDER, exist_ok=True)
    result_filename = os.path.basename(input_file)  # Get only the filename
    result_path = os.path.join(RESULT_FOLDER, result_filename)
    os.system(f'cp {input_file} {result_path}')
    print(f'cp {input_file} {result_path}')
    return result_filename  # Return just the filename

@app.route('/download/<filename>', methods=['GET'])
def download_file(filename):
    try:
        filename = secure_filename(filename)  # Sanitize the file name
        return send_from_directory(RESULT_FOLDER, filename, as_attachment=True)
    except FileNotFoundError:
        return jsonify({'error': 'File not found'}), 404

@app.route('/preview/<filename>', methods=['GET'])
def preview_file(filename):
    result_filepath = os.path.join(app.config['RESULT_FOLDER'], filename)
    if not os.path.exists(result_filepath):
        return jsonify({'error': 'File not found'}), 404

    with open(result_filepath, 'r') as result_file:
        content = result_file.read()

    return jsonify({'content': content})



if __name__ == '__main__':
    app.run(debug=True)