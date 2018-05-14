/*

Tailwind - The Utility-First CSS Framework

A project by Adam Wathan (@adamwathan), Jonathan Reinink (@reinink),
David Hemphill (@davidhemphill) and Steve Schoger (@steveschoger).

Welcome to the Tailwind config file. This is where you can customize
Tailwind specifically for your project. Don't be intimidated by the
length of this file. It's really just a big JavaScript object and
we've done our very best to explain each section.

View the full documentation at https://tailwindcss.com.


|-------------------------------------------------------------------------------
| The default config
|-------------------------------------------------------------------------------
|
| This variable contains the default Tailwind config. You don't have
| to use it, but it can sometimes be helpful to have available. For
| example, you may choose to merge your custom configuration
| values with some of the Tailwind defaults.
|
*/

let defaultConfig = require('tailwindcss/defaultConfig')()

let mycolors = {
  'bg-primary': 'transparent',
  'bg-secondary': 'var(--color-bg-secondary)',
  'bg-default': 'var(--color-bg-default)',
  'bg-inverse': '#3d4852',

  'text-primary': '#3d4852',
  'text-secondary': 'var(--color-text-secondary)',
  'text-default': 'var(--color-text-default)',
  'text-default-soft': 'var(--color-text-default-soft)',
  'text-inverse': '#f8fafc',
  'text-inverse-soft': 'var(--color-text-inverse-soft)',
}

module.exports = {
  mycolors: mycolors,

  fonts: {
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

  plugins: [
    require('tailwindcss/plugins/container')({
      center: true,
      // padding: '1rem',
    }),
  ],
}
