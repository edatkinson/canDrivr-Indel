

from app import temp_model
from werkzeug.utils import secure_filename
filename = 'Einstein’s Theories.md'
filename = secure_filename(filename)
print(filename)



