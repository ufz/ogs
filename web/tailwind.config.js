let primaryColor = "#104EB2";
let accentColor = "#f6993f";
let textColor = "#333333";
let textInverseColor = "#f8fafc"

let mycolors = {
  'bg-primary': 'transparent',
  'bg-secondary': 'var(--color-bg-secondary)',
  'bg-default': 'var(--color-bg-default)',
  'bg-inverse': primaryColor,

  'text-primary': primaryColor,
  'text-accent': accentColor,
  'text-default': textColor,
  'text-default-soft': 'var(--color-text-default-soft)',
  'text-inverse': textInverseColor,
  'text-inverse-soft': 'var(--color-text-inverse-soft)',
}

module.exports = {
  purge: [
    './layouts/**/*.html',
  ],
  theme: {
    mycolors: mycolors,
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
