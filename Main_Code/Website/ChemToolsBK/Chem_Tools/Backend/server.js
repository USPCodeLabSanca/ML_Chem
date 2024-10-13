const express = require('express');
const path = require('path');

const app = express();
const port = 3000;

// Serve static files from the React app build directory
app.use(express.static(path.join(__dirname, '../Frontend/dist')));

// Handle any requests that don't match the above
app.get('*', (req, res) => {
  res.sendFile(path.join(__dirname, '../Frontend/dist', 'index.html'));
});

app.listen(port, () => {
  console.log(`Server running at http://localhost:${port}/`);
});
