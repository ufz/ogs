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
        'brand': {
          50: '#F3F6FB',
          100: '#E7EDF7',
          200: '#C3D3EC',
          300: '#9FB8E0',
          400: '#5883C9',
          500: '#104EB2',
          600: '#0E46A0',
          700: '#0A2F6B',
          800: '#072350',
          900: '#051735',
        },
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
    typography: (theme) => ({
      default: {
        css: {
          a: {
            color: theme('colors.brand.500'),
          }
        },
      },
    }),
  },
  plugins: [
    require('@tailwindcss/typography'),
  ],
}
