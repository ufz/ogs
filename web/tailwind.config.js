module.exports = {
  future: {
    removeDeprecatedGapUtilities: true,
  },
  purge: [
    './layouts/**/*.html',
  ],
  theme: {
    extend: {
      fontSize: {
        '8xl': '6rem',
      },
      colors: {
        'brand': '#104EB2',
      },
      spacing: {
        '72': '18rem',
      }
    },
    fontFamily: {
      'sans': [
        'Open Sans',
        'system-ui',
        'BlinkMacSystemFont',
        '-apple-system',
        'Segoe UI',
        'Oxygen',
        'Ubuntu',
        'Cantarell',
        'Fira Sans',
        'Droid Sans',
        'Helvetica Neue',
        'sans-serif',
      ],
      'serif': [
        'Roboto',
        'Constantia',
        'Lucida Bright',
        'Lucidabright',
        'Lucida Serif',
        'Lucida',
        'DejaVu Serif',
        'Bitstream Vera Serif',
        'Liberation Serif',
        'Georgia',
        'serif',
      ],
      'serif-caps': [
        'serif'
      ],
      'mono': [
        'Menlo',
        'Monaco',
        'Consolas',
        'Liberation Mono',
        'Courier New',
        'monospace',
      ]
    },
  },
  plugins: [
    require('@tailwindcss/typography'),
  ],
}
