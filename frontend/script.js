
document.querySelector('form').addEventListener('submit', async (event) => {
  event.preventDefault();
  const formData = new FormData();
  const fileInput = document.getElementById('file');
  formData.append('file', fileInput.files[0]);

  const response = await fetch('/upload', {
      method: 'POST',
      body: formData,
  });

  const result = await response.json();

  if (result.result_file) {
      const downloadButton = document.getElementById('downloadButton');
      const previewDiv = document.getElementById('preview');
      
      // Enable the download button
      downloadButton.disabled = false;
      downloadButton.addEventListener('click', () => {
          const encodedFilename = encodeURIComponent(result.result_file);
          window.location.href = `/download/${encodedFilename}`;
      });

      // Fetch and display preview
      const previewResponse = await fetch(`/preview/${encodeURIComponent(result.result_file)}`);
      const previewData = await previewResponse.json();

      if (previewData.content) {
          previewDiv.innerText = previewData.content;
          previewDiv.style.display = 'block'; // Make sure it's visible
      } else {
          previewDiv.innerText = 'Error loading preview.';
          previewDiv.style.display = 'block';
      }
  } else {
      console.error('Error:', result.error);
  }
});


function noDownload() {
    const downloadButton = document.getElementById('downloadButton');
    downloadButton.disabled = true; // Disable the button
    downloadButton.onclick = () => alert('No file to download.');
    downloadButton.innerText = 'No file to download';
    downloadButton.title = 'Please upload a file first';
}

function yesDownload() {
    const downloadButton = document.getElementById('downloadButton');
    downloadButton.disabled = false; // Enable the button
    downloadButton.onclick = handleDownloadClick; // Attach the correct click handler
    downloadButton.innerText = 'Download';
    downloadButton.title = '';
}

function handleDownloadClick() {
    alert('File Downloaded!');
    // Add actual download logic here if needed
}

function isFileSelected() {
    const fileInput = document.getElementById('file');
    return fileInput.files && fileInput.files.length > 0; // Check if a file is selected
}

// Attach an event listener to the file input to check file selection
document.getElementById('file').addEventListener('change', () => {
    if (isFileSelected()) {
        yesDownload();
    } else {
        noDownload();
    }
});

// Initial setup: Ensure the download button is disabled when the page loads
document.addEventListener('DOMContentLoaded', () => {
    noDownload();
});

window.addEventListener('scroll', () => {
    const scrollPosition = window.scrollY; // Get the current vertical scroll position
    const documentHeight = document.documentElement.scrollHeight - window.innerHeight;

    // Calculate the scroll percentage (0 to 1)
    const scrollPercentage = scrollPosition / documentHeight;

    // Calculate the RGB values for the blue gradient
    const red = Math.floor(220 - (220 - 70) * scrollPercentage); // Transition from 173 to 70
    const green = Math.floor(220 - (220 - 130) * scrollPercentage); // Transition from 216 to 130
    const blue = Math.floor(230 - (230 - 200) * scrollPercentage); // Transition from 230 to 180

    // Apply the color to the body
    document.body.style.backgroundColor = `rgb(${red}, ${green}, ${blue})`;
});
